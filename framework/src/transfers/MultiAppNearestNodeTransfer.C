//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiAppNearestNodeTransfer.h"

// MOOSE includes
#include "DisplacedProblem.h"
#include "FEProblem.h"
#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseVariableFE.h"

#include "libmesh/system.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_algebra.h"

registerMooseObject("MooseApp", MultiAppNearestNodeTransfer);

template <>
InputParameters
validParams<MultiAppNearestNodeTransfer>()
{
  InputParameters params = validParams<MultiAppTransfer>();

  params.addRequiredParam<AuxVariableName>(
      "variable", "The auxiliary variable to store the transferred values in.");
  params.addRequiredParam<VariableName>("source_variable", "The variable to transfer from.");
  params.addParam<BoundaryName>(
      "source_boundary",
      "The boundary we are transferring from (if not specified, whole domain is used).");
  params.addParam<BoundaryName>(
      "target_boundary",
      "The boundary we are transferring to (if not specified, whole domain is used).");
  params.addParam<bool>("fixed_meshes",
                        false,
                        "Set to true when the meshes are not changing (ie, "
                        "no movement or adaptivity).  This will cache "
                        "nearest node neighbors to greatly speed up the "
                        "transfer.");

  return params;
}

MultiAppNearestNodeTransfer::MultiAppNearestNodeTransfer(const InputParameters & parameters)
  : MultiAppTransfer(parameters),
    _to_var_name(getParam<AuxVariableName>("variable")),
    _from_var_name(getParam<VariableName>("source_variable")),
    _fixed_meshes(getParam<bool>("fixed_meshes")),
    _node_map(declareRestartableData<std::map<dof_id_type, Node *>>("node_map")),
    _distance_map(declareRestartableData<std::map<dof_id_type, Real>>("distance_map")),
    _neighbors_cached(declareRestartableData<bool>("neighbors_cached", false)),
    _cached_froms(declareRestartableData<std::vector<std::vector<unsigned int>>>("cached_froms")),
    _cached_dof_ids(
        declareRestartableData<std::vector<std::vector<dof_id_type>>>("cached_dof_ids")),
    _cached_from_inds(
        declareRestartableData<std::map<dof_id_type, unsigned int>>("cached_from_ids")),
    _cached_qp_inds(declareRestartableData<std::map<dof_id_type, unsigned int>>("cached_qp_inds"))
{
}

void
MultiAppNearestNodeTransfer::initialSetup()
{
  if (_direction == TO_MULTIAPP)
    variableIntegrityCheck(_to_var_name);
  else
    variableIntegrityCheck(_from_var_name);
}

void
MultiAppNearestNodeTransfer::execute()
{
  _console << "Beginning NearestNodeTransfer " << name() << std::endl;

  getAppInfo();

  // Get the bounding boxes for the "from" domains.
  std::vector<BoundingBox> bboxes = getFromBoundingBoxes();

  // Debugging:
  for (const auto & bbox : bboxes)
    std::cout << "bbox = (" << bbox.first << ", " << bbox.second << ")" << std::endl;

  // Figure out how many "from" domains each processor owns.
  std::vector<unsigned int> froms_per_proc = getFromsPerProc();

  // Debugging:
  std::cout << "froms_per_proc = ";
  for (const auto & val : froms_per_proc)
    std::cout << val << " ";
  std::cout << std::endl;

  ////////////////////
  // For every point in the local "to" domain, figure out which "from" domains
  // might contain it's nearest neighbor, and send that point to the processors
  // that own those "from" domains.
  //
  // How do we know which "from" domains might contain the nearest neighbor, you
  // ask?  Well, consider two "from" domains, A and B.  If every point in A is
  // closer than every point in B, then we know that B cannot possibly contain
  // the nearest neighbor.  Hence, we'll only check A for the nearest neighbor.
  // We'll use the functions bboxMaxDistance and bboxMinDistance to figure out
  // if every point in A is closer than every point in B.
  ////////////////////

  // outgoing_qps = nodes/centroids we'll send to other processors.
  std::vector<std::vector<Point>> outgoing_qps(n_processors());
  // When we get results back, node_index_map will tell us which results go with
  // which points
  std::vector<std::map<std::pair<unsigned int, unsigned int>, unsigned int>> node_index_map(
      n_processors());

  if (!_neighbors_cached)
  {
    for (unsigned int i_to = 0; i_to < _to_problems.size(); i_to++)
    {
      System * to_sys = find_sys(*_to_es[i_to], _to_var_name);
      unsigned int sys_num = to_sys->number();
      unsigned int var_num = to_sys->variable_number(_to_var_name);
      MeshBase * to_mesh = &_to_meshes[i_to]->getMesh();
      bool is_nodal = to_sys->variable_type(var_num).family == LAGRANGE;

      if (is_nodal)
      {
        std::vector<Node *> target_local_nodes;

        if (isParamValid("target_boundary"))
        {
          BoundaryID target_bnd_id =
              _to_meshes[i_to]->getBoundaryID(getParam<BoundaryName>("target_boundary"));

          ConstBndNodeRange & bnd_nodes = *(_to_meshes[i_to])->getBoundaryNodeRange();
          for (const auto & bnode : bnd_nodes)
            if (bnode->_bnd_id == target_bnd_id && bnode->_node->processor_id() == processor_id())
              target_local_nodes.push_back(bnode->_node);
        }
        else
        {
          target_local_nodes.resize(to_mesh->n_local_nodes());
          unsigned int i = 0;
          for (auto & node : to_mesh->local_node_ptr_range())
            target_local_nodes[i++] = node;
        }

        std::cout << "Finished building list of target_local_nodes: ";
        for (const auto & node : target_local_nodes)
          std::cout << node->id() << " ";
        std::cout << std::endl;

        for (const auto & node : target_local_nodes)
        {
          // Newline at the start of each node iteration
          std::cout << std::endl;

          // Skip this node if the variable has no dofs at it.
          if (node->n_dofs(sys_num, var_num) < 1)
            continue;

          // Find which bboxes might have the nearest node to this point.
          Real nearest_max_distance = std::numeric_limits<Real>::max();
          for (const auto & bbox : bboxes)
          {
            Real distance = bboxMaxDistance(*node, bbox);
            if (distance < nearest_max_distance)
              nearest_max_distance = distance;
          }

          std::cout << "The nearest_max_distance (over all bboxes) for Node " << node->id()
                    << ", " << static_cast<Point&>(*node) << " = " << nearest_max_distance << std::endl;

          unsigned int from0 = 0;
          for (processor_id_type i_proc = 0; i_proc < n_processors();
               from0 += froms_per_proc[i_proc], i_proc++)
          {
            bool qp_found = false;
            for (unsigned int i_from = from0; i_from < from0 + froms_per_proc[i_proc] && !qp_found;
                 i_from++)
            {

              Real distance = bboxMinDistance(*node, bboxes[i_from]);
              std::cout << "Testing bboxMinDistance = " << distance << " to bbox i_from = " << i_from << std::endl;

              if (distance <= nearest_max_distance || bboxes[i_from].contains_point(*node))
              {
                // Debugging
                std::cout << "Setting node_index_map info for node " << node->id() << std::endl;

                // Which condition above actually got us here?
                if (distance <= nearest_max_distance)
                  std::cout << "distance = " << distance << " less than or equal to " << nearest_max_distance << std::endl;

                if (bboxes[i_from].contains_point(*node))
                  std::cout << "bboxes[" << i_from << "] contains point " << static_cast<Point&>(*node) << std::endl;
                else
                  std::cout << "bboxes[" << i_from << "] did not contain point " << static_cast<Point&>(*node) << std::endl;

                std::pair<unsigned int, unsigned int> key(i_to, node->id());
                std::cout << "i_proc = " << i_proc << std::endl;
                std::cout << "key = (i_to, node->id()) = " << i_to << ", " << node->id() << std::endl;
                std::cout << "outgoing_qps[i_proc].size() = " << outgoing_qps[i_proc].size() << std::endl;

                node_index_map[i_proc][key] = outgoing_qps[i_proc].size();
                outgoing_qps[i_proc].push_back(*node + _to_positions[i_to]);
                qp_found = true;
              }

              // Debugging: if (!qp_found), this means something went
              // wrong. It can happen if
              // bboxMinDistance == nearest_max_distance
              // and the point is not contained in the bbox. This should not happen in principle, since
              // for a given BBox A, bboxMinDistance(p,A) <= bboxMaxDistance(p,A). So perhaps changing
              // the check to <= would be sufficient...
              if (!qp_found)
                {
                  std::cout << "Node " << node->id() << ", " << static_cast<Point&>(*node) << " not found closest to _any_ BBox." << std::endl;
                  // Based on the logic above, we should have found at
                  // least one bounding box for every
                  // target_local_node, since at least one min
                  // distance will be less than or equal to the
                  // nearest_max_distance, at the very least this must
                  // happen for the bbox which sets the
                  // nearest_max_distance.
                  mooseError("BoundingBox found for node ", node->id(), " not found.");
                }
            }
          }
        }
      }
      else // Elemental
      {
        for (auto & elem : as_range(to_mesh->local_elements_begin(), to_mesh->local_elements_end()))
        {
          Point centroid = elem->centroid();

          // Skip this element if the variable has no dofs at it.
          if (elem->n_dofs(sys_num, var_num) < 1)
            continue;

          // Find which bboxes might have the nearest node to this point.
          Real nearest_max_distance = std::numeric_limits<Real>::max();
          for (const auto & bbox : bboxes)
          {
            Real distance = bboxMaxDistance(centroid, bbox);
            if (distance < nearest_max_distance)
              nearest_max_distance = distance;
          }

          unsigned int from0 = 0;
          for (processor_id_type i_proc = 0; i_proc < n_processors();
               from0 += froms_per_proc[i_proc], i_proc++)
          {
            bool qp_found = false;
            for (unsigned int i_from = from0; i_from < from0 + froms_per_proc[i_proc] && !qp_found;
                 i_from++)
            {
              Real distance = bboxMinDistance(centroid, bboxes[i_from]);
              if (distance < nearest_max_distance || bboxes[i_from].contains_point(centroid))
              {
                std::pair<unsigned int, unsigned int> key(i_to, elem->id());
                node_index_map[i_proc][key] = outgoing_qps[i_proc].size();
                outgoing_qps[i_proc].push_back(centroid + _to_positions[i_to]);
                qp_found = true;
              }
            }
          }
        }
      }
    }
    // We should be done building the node_index_map and outgoing_qps data structures.
    // node_index_map = vector<map<pair<unsigned int, unsigned int>, unsigned int>>
    unsigned int pid=0;
    std::cout << "Printing entries of node_index_map:" << std::endl;
    for (const auto & m : node_index_map)
    {
      std::cout << "Map entries for processor " << pid++ << std::endl;
      for (const auto & pr : m)
      {
        const auto & key = pr.first;
        // std::cout << "Key: (send to proc = " << key.first << ", node id = " << key.second << ")" << std::endl;
        // std::cout << "Index of point: " << pr.second << std::endl;
        std::cout << "[" << pid << "][(" << key.first << ", " << key.second << ")] = " << pr.second << std::endl;
      }
    }
  }

  ////////////////////
  // Send local node/centroid positions off to the other processors and take
  // care of points sent to this processor.  We'll need to check the points
  // against all of the "from" domains that this processor owns.  For each
  // point, we'll find the nearest node, then we'll send the value at that node
  // and the distance between the node and the point back to the processor that
  // requested that point.
  ////////////////////

  std::vector<std::vector<Real>> incoming_evals(n_processors());
  std::vector<Parallel::Request> send_qps(n_processors());
  std::vector<Parallel::Request> send_evals(n_processors());

  // Create these here so that they live the entire life of this function
  // and are NOT reused per processor.
  std::vector<std::vector<Real>> processor_outgoing_evals(n_processors());

  if (!_neighbors_cached)
  {
    for (processor_id_type i_proc = 0; i_proc < n_processors(); i_proc++)
    {
      if (i_proc == processor_id())
        continue;
      _communicator.send(i_proc, outgoing_qps[i_proc], send_qps[i_proc]);
    }

    // Build an array of pointers to all of this processor's local nodes.  We
    // need to do this to avoid the expense of using LibMesh iterators.  This
    // step also takes care of limiting the search to boundary nodes, if
    // applicable.
    std::vector<std::vector<Node *>> local_nodes(froms_per_proc[processor_id()]);
    for (unsigned int i = 0; i < froms_per_proc[processor_id()]; i++)
    {
      getLocalNodes(_from_meshes[i], local_nodes[i]);
    }

    if (_fixed_meshes)
    {
      _cached_froms.resize(n_processors());
      _cached_dof_ids.resize(n_processors());
    }

    for (processor_id_type i_proc = 0; i_proc < n_processors(); i_proc++)
    {
      std::vector<Point> incoming_qps;
      if (i_proc == processor_id())
        incoming_qps = outgoing_qps[i_proc];
      else
        _communicator.receive(i_proc, incoming_qps);

      if (_fixed_meshes)
      {
        _cached_froms[i_proc].resize(incoming_qps.size());
        _cached_dof_ids[i_proc].resize(incoming_qps.size());
      }

      std::vector<Real> & outgoing_evals = processor_outgoing_evals[i_proc];
      outgoing_evals.resize(2 * incoming_qps.size());

      for (unsigned int qp = 0; qp < incoming_qps.size(); qp++)
      {
        Point qpt = incoming_qps[qp];
        outgoing_evals[2 * qp] = std::numeric_limits<Real>::max();
        for (unsigned int i_local_from = 0; i_local_from < froms_per_proc[processor_id()];
             i_local_from++)
        {
          MooseVariableFEBase & from_var =
              _from_problems[i_local_from]->getVariable(0,
                                                        _from_var_name,
                                                        Moose::VarKindType::VAR_ANY,
                                                        Moose::VarFieldType::VAR_FIELD_STANDARD);
          System & from_sys = from_var.sys().system();
          unsigned int from_sys_num = from_sys.number();
          unsigned int from_var_num = from_sys.variable_number(from_var.name());

          for (unsigned int i_node = 0; i_node < local_nodes[i_local_from].size(); i_node++)
          {
            Real current_distance =
                (qpt - *(local_nodes[i_local_from][i_node]) - _from_positions[i_local_from]).norm();
            if (current_distance < outgoing_evals[2 * qp])
            {
              // Assuming LAGRANGE!
              if (local_nodes[i_local_from][i_node]->n_dofs(from_sys_num, from_var_num) > 0)
              {
                dof_id_type from_dof =
                    local_nodes[i_local_from][i_node]->dof_number(from_sys_num, from_var_num, 0);

                outgoing_evals[2 * qp] = current_distance;
                outgoing_evals[2 * qp + 1] = (*from_sys.solution)(from_dof);

                if (_fixed_meshes)
                {
                  // Cache the nearest nodes.
                  _cached_froms[i_proc][qp] = i_local_from;
                  _cached_dof_ids[i_proc][qp] = from_dof;
                }
              }
            }
          }
        }
      }

      if (i_proc == processor_id())
        incoming_evals[i_proc] = outgoing_evals;
      else
        _communicator.send(i_proc, outgoing_evals, send_evals[i_proc]);
    }

    // Every proc should now have filled in the "processor_outgoing_evals" data structure (?)
    std::cout << "processor_outgoing_evals = " << std::endl;
    for (unsigned int pid=0; pid<processor_outgoing_evals.size(); ++pid) // const auto & vec : processor_outgoing_evals)
      for (unsigned int i=0; i<processor_outgoing_evals[pid].size(); ++i) // const auto & val : vec)
        std::cout << "[" << pid << "][" << i << "] = " << processor_outgoing_evals[pid][i] << std::endl;
  }

  else // We've cached the nearest nodes.
  {
    for (processor_id_type i_proc = 0; i_proc < n_processors(); i_proc++)
    {
      std::vector<Real> & outgoing_evals = processor_outgoing_evals[i_proc];
      outgoing_evals.resize(_cached_froms[i_proc].size());

      for (unsigned int qp = 0; qp < outgoing_evals.size(); qp++)
      {
        MooseVariableFEBase & from_var = _from_problems[_cached_froms[i_proc][qp]]->getVariable(
            0,
            _from_var_name,
            Moose::VarKindType::VAR_ANY,
            Moose::VarFieldType::VAR_FIELD_STANDARD);
        System & from_sys = from_var.sys().system();
        dof_id_type from_dof = _cached_dof_ids[i_proc][qp];
        // outgoing_evals[qp] = (*from_sys.solution)(_cached_dof_ids[i_proc][qp]);
        outgoing_evals[qp] = (*from_sys.solution)(from_dof);
      }

      if (i_proc == processor_id())
        incoming_evals[i_proc] = outgoing_evals;
      else
        _communicator.send(i_proc, outgoing_evals, send_evals[i_proc]);
    }
  }

  ////////////////////
  // Gather all of the evaluations, find the nearest one for each node/element,
  // and apply the values.
  ////////////////////

  for (processor_id_type i_proc = 0; i_proc < n_processors(); i_proc++)
  {
    if (i_proc == processor_id())
      continue;

    _communicator.receive(i_proc, incoming_evals[i_proc]);
  }

  for (unsigned int i_to = 0; i_to < _to_problems.size(); i_to++)
  {
    // Loop over the master nodes and set the value of the variable
    System * to_sys = find_sys(*_to_es[i_to], _to_var_name);

    unsigned int sys_num = to_sys->number();
    unsigned int var_num = to_sys->variable_number(_to_var_name);

    NumericVector<Real> * solution = nullptr;
    switch (_direction)
    {
      case TO_MULTIAPP:
        solution = &getTransferVector(i_to, _to_var_name);
        break;
      case FROM_MULTIAPP:
        solution = to_sys->solution.get();
        break;
      default:
        mooseError("Unknown direction");
    }

    MeshBase * to_mesh = &_to_meshes[i_to]->getMesh();

    bool is_nodal = to_sys->variable_type(var_num).family == LAGRANGE;

    if (is_nodal)
    {
      std::vector<Node *> target_local_nodes;

      if (isParamValid("target_boundary"))
      {
        BoundaryID target_bnd_id =
            _to_meshes[i_to]->getBoundaryID(getParam<BoundaryName>("target_boundary"));

        ConstBndNodeRange & bnd_nodes = *(_to_meshes[i_to])->getBoundaryNodeRange();
        for (const auto & bnode : bnd_nodes)
          if (bnode->_bnd_id == target_bnd_id && bnode->_node->processor_id() == processor_id())
            target_local_nodes.push_back(bnode->_node);
      }
      else
      {
        target_local_nodes.resize(to_mesh->n_local_nodes());
        unsigned int i = 0;
        for (auto & node : to_mesh->local_node_ptr_range())
          target_local_nodes[i++] = node;
      }

      for (const auto & node : target_local_nodes)
      {
        // Skip this node if the variable has no dofs at it.
        if (node->n_dofs(sys_num, var_num) < 1)
          continue;

        Real best_val = 0;
        if (!_neighbors_cached)
        {
          Real min_dist = std::numeric_limits<Real>::max();
          for (unsigned int i_from = 0; i_from < incoming_evals.size(); i_from++)
          {
            std::pair<unsigned int, unsigned int> key(i_to, node->id());
            if (node_index_map[i_from].find(key) == node_index_map[i_from].end())
              continue;
            unsigned int qp_ind = node_index_map[i_from][key];
            if (incoming_evals[i_from][2 * qp_ind] >= min_dist)
              continue;
            min_dist = incoming_evals[i_from][2 * qp_ind];
            best_val = incoming_evals[i_from][2 * qp_ind + 1];

            if (_fixed_meshes)
            {
              // Cache these indices.
              _cached_from_inds[node->id()] = i_from;
              _cached_qp_inds[node->id()] = qp_ind;
            }
          }
        }

        else
        {
          best_val = incoming_evals[_cached_from_inds[node->id()]][_cached_qp_inds[node->id()]];
        }

        dof_id_type dof = node->dof_number(sys_num, var_num, 0);
        solution->set(dof, best_val);
      }
    }
    else // Elemental
    {
      for (auto & elem : as_range(to_mesh->local_elements_begin(), to_mesh->local_elements_end()))
      {
        // Skip this element if the variable has no dofs at it.
        if (elem->n_dofs(sys_num, var_num) < 1)
          continue;

        Real best_val = 0;
        if (!_neighbors_cached)
        {
          Real min_dist = std::numeric_limits<Real>::max();
          for (unsigned int i_from = 0; i_from < incoming_evals.size(); i_from++)
          {
            std::pair<unsigned int, unsigned int> key(i_to, elem->id());
            if (node_index_map[i_from].find(key) == node_index_map[i_from].end())
              continue;
            unsigned int qp_ind = node_index_map[i_from][key];
            if (incoming_evals[i_from][2 * qp_ind] >= min_dist)
              continue;
            min_dist = incoming_evals[i_from][2 * qp_ind];
            best_val = incoming_evals[i_from][2 * qp_ind + 1];

            if (_fixed_meshes)
            {
              // Cache these indices.
              _cached_from_inds[elem->id()] = i_from;
              _cached_qp_inds[elem->id()] = qp_ind;
            }
          }
        }

        else
        {
          best_val = incoming_evals[_cached_from_inds[elem->id()]][_cached_qp_inds[elem->id()]];
        }

        dof_id_type dof = elem->dof_number(sys_num, var_num, 0);
        solution->set(dof, best_val);
      }
    }
    solution->close();
    to_sys->update();
  }

  if (_fixed_meshes)
    _neighbors_cached = true;

  // Make sure all our sends succeeded.
  for (processor_id_type i_proc = 0; i_proc < n_processors(); i_proc++)
  {
    if (i_proc == processor_id())
      continue;
    send_qps[i_proc].wait();
    send_evals[i_proc].wait();
  }

  _console << "Finished NearestNodeTransfer " << name() << std::endl;
}

Node *
MultiAppNearestNodeTransfer::getNearestNode(const Point & p,
                                            Real & distance,
                                            MooseMesh * mesh,
                                            bool local)
{
  distance = std::numeric_limits<Real>::max();
  Node * nearest = NULL;

  if (isParamValid("source_boundary"))
  {
    BoundaryID src_bnd_id = mesh->getBoundaryID(getParam<BoundaryName>("source_boundary"));

    ConstBndNodeRange & bnd_nodes = *mesh->getBoundaryNodeRange();
    for (const auto & bnode : bnd_nodes)
    {
      if (bnode->_bnd_id == src_bnd_id)
      {
        Node * node = bnode->_node;
        Real current_distance = (p - *node).norm();

        if (current_distance < distance)
        {
          distance = current_distance;
          nearest = node;
        }
      }
    }
  }
  else
  {
    for (auto & node : as_range(local ? mesh->localNodesBegin() : mesh->getMesh().nodes_begin(),
                                local ? mesh->localNodesEnd() : mesh->getMesh().nodes_end()))
    {
      Real current_distance = (p - *node).norm();

      if (current_distance < distance)
      {
        distance = current_distance;
        nearest = node;
      }
    }
  }

  return nearest;
}

Real
MultiAppNearestNodeTransfer::bboxMaxDistance(Point p, BoundingBox bbox)
{
  std::vector<Point> source_points = {bbox.first, bbox.second};

  // std::cout << "Computing max distance between Point " << p << " and BoundingBox (" << bbox.first << ", " << bbox.second << ")" << std::endl;

  std::vector<Point> all_points(8);
  for (unsigned int x = 0; x < 2; x++)
    for (unsigned int y = 0; y < 2; y++)
      for (unsigned int z = 0; z < 2; z++)
        all_points[x + 2 * y + 4 * z] =
            Point(source_points[x](0), source_points[y](1), source_points[z](2));

  Real max_distance = 0.;

  // Compute the signed distance from p to bbox.
  // Real signed_distance = bbox.signed_distance(p);
  // std::cout << "The signed distance between Point " << p << " and BoundingBox (" << bbox.first << ", " << bbox.second << ") = " << signed_distance << std::endl;

  for (unsigned int i = 0; i < 8; i++)
  {
    Real distance = (p - all_points[i]).norm();
    if (distance > max_distance)
      max_distance = distance;
  }

  // std::cout << "Returning max_distance=" << max_distance << std::endl;

  return max_distance;
}

Real
MultiAppNearestNodeTransfer::bboxMinDistance(Point p, BoundingBox bbox)
{
  std::vector<Point> source_points = {bbox.first, bbox.second};

  // std::cout << "Computing min distance between Point " << p << " and BoundingBox (" << bbox.first << ", " << bbox.second << ")" << std::endl;

  std::vector<Point> all_points(8);
  for (unsigned int x = 0; x < 2; x++)
    for (unsigned int y = 0; y < 2; y++)
      for (unsigned int z = 0; z < 2; z++)
        all_points[x + 2 * y + 4 * z] =
            Point(source_points[x](0), source_points[y](1), source_points[z](2));

  Real min_distance = std::numeric_limits<Real>::max();

  // Compute the signed distance from p to bbox.
  // Real signed_distance = bbox.signed_distance(p);
  // std::cout << "The signed distance between Point " << p << " and BoundingBox (" << bbox.first << ", " << bbox.second << ") = " << signed_distance << std::endl;

  for (unsigned int i = 0; i < 8; i++)
  {
    Real distance = (p - all_points[i]).norm();
    if (distance < min_distance)
      min_distance = distance;
  }

  // std::cout << "Returning min_distance=" << min_distance << std::endl;

  return min_distance;
}

void
MultiAppNearestNodeTransfer::getLocalNodes(MooseMesh * mesh, std::vector<Node *> & local_nodes)
{
  if (isParamValid("source_boundary"))
  {
    BoundaryID src_bnd_id = mesh->getBoundaryID(getParam<BoundaryName>("source_boundary"));

    ConstBndNodeRange & bnd_nodes = *mesh->getBoundaryNodeRange();
    for (const auto & bnode : bnd_nodes)
      if (bnode->_bnd_id == src_bnd_id && bnode->_node->processor_id() == processor_id())
        local_nodes.push_back(bnode->_node);
  }
  else
  {
    local_nodes.resize(mesh->getMesh().n_local_nodes());
    unsigned int i = 0;
    for (auto & node : as_range(mesh->localNodesBegin(), mesh->localNodesEnd()))
      local_nodes[i++] = node;
  }
}

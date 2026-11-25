#pragma once

#include "RelationshipManager.h"
#include <unordered_map>

/**
 * Ghosting and sparsity for non-neighbour connections (NNCs) between cells that do not share a
 * mesh face, e.g. across faults with throw in a corner-point grid. The pairs are read from a CSV
 * file with columns elem_a, elem_b, trans (element ids of the mesh as read, mesh renumbering off).
 */
class NNCRelationshipManager : public RelationshipManager
{
public:
  static InputParameters validParams();
  NNCRelationshipManager(const InputParameters & params);
  NNCRelationshipManager(const NNCRelationshipManager & other);

  void operator()(const MeshBase::const_element_iterator & range_begin,
                  const MeshBase::const_element_iterator & range_end,
                  processor_id_type p,
                  map_type & coupled_elements) override;

  std::unique_ptr<libMesh::GhostingFunctor> clone() const override;
  std::string getInfo() const override;
  bool operator>=(const RelationshipManager & other) const override;

  /// Reads an NNC file into a symmetric map elem -> [(partner, trans)]
  static void readNNCFile(const std::string & file,
                          std::unordered_map<dof_id_type, std::vector<std::pair<dof_id_type, Real>>> & map);

protected:
  const std::string _nnc_file;
  std::unordered_map<dof_id_type, std::vector<std::pair<dof_id_type, Real>>> _nnc;
};

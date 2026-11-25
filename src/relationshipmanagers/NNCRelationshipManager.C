#include "NNCRelationshipManager.h"
#include "MooseApp.h"
#include "MooseUtils.h"
#include "DelimitedFileReader.h"

#include "libmesh/elem.h"
#include "libmesh/mesh_base.h"

registerMooseObject("PorousFlowApp", NNCRelationshipManager);

InputParameters
NNCRelationshipManager::validParams()
{
  InputParameters params = RelationshipManager::validParams();
  params.addRequiredParam<FileName>("nnc_file", "CSV file: elem_a, elem_b, trans");
  params.set<bool>("attach_geometric_early") = false;
  return params;
}

NNCRelationshipManager::NNCRelationshipManager(const InputParameters & params)
  : RelationshipManager(params), _nnc_file(getParam<FileName>("nnc_file"))
{
  readNNCFile(_nnc_file, _nnc);
}

NNCRelationshipManager::NNCRelationshipManager(const NNCRelationshipManager & other)
  : RelationshipManager(other), _nnc_file(other._nnc_file), _nnc(other._nnc)
{
}

void
NNCRelationshipManager::readNNCFile(
    const std::string & file,
    std::unordered_map<dof_id_type, std::vector<std::pair<dof_id_type, Real>>> & map)
{
  MooseUtils::DelimitedFileReader reader(file);
  reader.setFormatFlag(MooseUtils::DelimitedFileReader::FormatFlag::COLUMNS);
  reader.read();
  const auto & a = reader.getData("elem_a");
  const auto & b = reader.getData("elem_b");
  const auto & t = reader.getData("trans");
  for (const auto i : index_range(a))
  {
    const auto ea = static_cast<dof_id_type>(std::lround(a[i]));
    const auto eb = static_cast<dof_id_type>(std::lround(b[i]));
    map[ea].emplace_back(eb, t[i]);
    map[eb].emplace_back(ea, t[i]);
  }
}

std::string
NNCRelationshipManager::getInfo() const
{
  return "NNCRelationshipManager (" + _nnc_file + ")";
}

void
NNCRelationshipManager::operator()(const MeshBase::const_element_iterator & range_begin,
                                   const MeshBase::const_element_iterator & range_end,
                                   const processor_id_type p,
                                   map_type & coupled_elements)
{
  static const libMesh::CouplingMatrix * const null_mat = nullptr;
  for (const Elem * const elem : as_range(range_begin, range_end))
  {
    const auto it = _nnc.find(elem->id());
    if (it == _nnc.end())
      continue;
    for (const auto & pr : it->second)
    {
      const Elem * other = _mesh->query_elem_ptr(pr.first);
      if (other && (other->processor_id() != p || _rm_type != Moose::RelationshipManagerType::ALGEBRAIC))
        coupled_elements.emplace(other, null_mat);
    }
  }
}

bool
NNCRelationshipManager::operator>=(const RelationshipManager & other) const
{
  if (auto * o = dynamic_cast<const NNCRelationshipManager *>(&other))
    return _nnc_file == o->_nnc_file && baseGreaterEqual(other);
  return false;
}

std::unique_ptr<libMesh::GhostingFunctor>
NNCRelationshipManager::clone() const
{
  return _app.getFactory().copyConstruct(*this);
}

half_aq_thickness = 25

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    xmin = 0.1
    dx = '1 2 4 43'
    ix = '10 10 10 10'
    dy = '30 10 40 10 30'
    iy = '10 10 10 10 10'
  []
  [shift]
    type = TransformGenerator
    input = mesh
    transform = translate
    vector_value = '0.1 -60 0'
  []
  [subdomains]
    type = ParsedSubdomainMeshGenerator
    input = shift
    combinatorial_geometry = 'y >= -${half_aq_thickness} & y <= ${half_aq_thickness}'
    block_id = 1
  []
  [rename]
    type = RenameBlockGenerator
    old_block = '0 1'
    new_block = 'caps aquifer'
    input = subdomains
  []
  coord_type = RZ
[]


# FEBio-preCICE adapter

**experimental** preCICE-adapter for the software tool [FEBio](https://github.com/febiosoftware/FEBio).

The adapter allows to exchange data on `FEMaterialPoints` with other solvers. The adapter is configured using a JSON file.

## Installing the plugin

We provide a x86 shared library as download option. This `.so` file requires preCICE, rttr and FEBio3.2 to be installed.

## Building the plugin

To build this plugin, get a checkout of the FEBio version you want to build with. Set the cmake variable `FEBio_LIB_DIR` to this dir, as we need the include files to compile this plugin.
Then run:

```
cmake
make
```

this produces the `libFEBioPreciceAdapter.so` in your cmake output folder. This is the file you need to reference in your `febio.xml`.

## Loading the plugin
Create a `febio.xml`

```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_config version="3.0">
        <default_linear_solver type="pardiso"></default_linear_solver>
        <import>$pathToThePluginFile/libFEBioPreciceAdapter.so</import>
</febio_config>
```

in your model folder.
Add the following sinppet to your `your_input_file.feb`  right below the Module section.

```xml
        <Code>
                <callback name="precice_callback"/>
        </Code>
```

Add a preCICE config file and create a `febio-config.json` with the options below
Than call febio with `febio3 -cnf ./febio.xml you_input_file.feb`

# configuration options
all attributes are contained in the `coupling_params` dict
| Key | Possible Values | Description |
| --- | --- | --- |
| `element_set_to_couple` |  | Name of the element set that should be coupled via preCICE. Can be found in the febio input file. |
| `participant_name`  |  | Name of the preCICE participant of FEBio |
| `config_file_name`  |  |  Path to the preCICE configuration file  |
| `read_mesh_name` |  | Name of the mesh to read from |
| `read_data_name` |  | dict, see below |
| `write_mesh_name` |  | name of the mesh to write to |
| `write_data_name` |  | dict, see below |

## write/read_data_name options
| Key | Possible Values | Description |
| --- | --- | --- |
| `type` | `function` `variable` | type of attribute for the FEMaterialPoint that is read or written. |
| `class_type` | `materialPoint` | This option sets what class type is accessed in FEBio, in future the possiblity of a `FEMaterial` or something like it could be added. |
| `febio_class_name` |  | Class name of the object that is being manipulated in FEBio e.g `FEFluidMaterialPoint` |
| `name` |  | Name of the attribute that is written/read `Vext` |
| `febio_type` | `vec3d` `vec2d` `double` `int` `float` `vector<dobule>` | data type of attribute being write/read |
| `mapping_name` |  | Name of data point in preCICE |
| `precice_type` | `scalar` `vector` | type of data point in precice |

## Example configuration

```json
{
   "coupling_params": {
      "element_set_to_couple": "Part1",
      "participant_name": "FEBio",
      "config_file_name": "../precice-config.xml",
      "read_mesh_name": "FEBioMesh",
      "read_data_name": [
         {"type": "function", "febio_class_name": "FESolutesMaterialPointTPM", "name": "set_concentrations", "mapping_name": "apap_ext'", "febio_type": "double", "precice_type": "scalar"},
         {"type": "function", "febio_class_name": "FESolutesMaterialPointTPM", "name": "set_concentrations_tangents", "mapping_name":  "d_vapap_ext__d_apap_ext", "febio_type": "double", "precice_type": "scalar"}
      ],
      "write_mesh_name": "FEBioMesh",
      "write_data_name": [
         {"type": "function", "class_type": "materialPoint", "febio_class_name": "FESolutesMaterialPointTPM", "name": "Vext", "mapping_name": "Vext", "febio_type": "double", "precice_type": "scalar"},
         {"type": "function", "class_type": "materialPoint", "febio_class_name": "FESolutesMaterialPointTPM", "name": "GetViFat", "mapping_name": "Vli_fat", "febio_type": "double", "precice_type": "scalar"},
         {"type": "variable", "class_type": "materialPoint", "febio_class_name": "FEBiphasicMaterialPointTPM", "name": "m_phi0", "mapping_name": "Vli_nofat", "febio_type": "double", "precice_type": "scalar"},
         {"type": "function", "class_type": "materialPoint", "febio_class_name": "FESolutesMaterialPointTPM", "name": "GetCyp2e1", "mapping_name": "cyp2e1", "febio_type": "double", "precice_type": "scalar"},
         {"type": "variable", "class_type": "materialPoint", "febio_class_name": "FESolutesMaterialPointTPM", "name": "m_ca", "mapping_name": "apap_ext", "febio_type": "vector<double>", "precice_type": "scalar"}
      ]
   }
}
```

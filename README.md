# FEBio-preCICE adapter

**experimental** preCICE-adapter for the software tool [FEBio](https://github.com/febiosoftware/FEBio).

The FEBio-preCICE adapter is a preCICE adapter for [FEBio](https://github.com/febiosoftware/FEBio). It allows to exchange data on `FEMaterialPoints`  via preCICE. What data and to which preCICE variable is configured via a JSON file (see example file `febio-config.json` below). For an easy setup use the provided `.so` file, or follow the build step to compile it for your platform. The adapter allows to exchange data on `FEMaterialPoints` with other solvers. The adapter is configured using a JSON file.

Design and development of this adapter is described in detail in the [Master thesis of Fritz Otlinghaus](http://dx.doi.org/10.18419/opus-12291).

## Required Dependencies

The following dependencies are required:

- [preCICE](https://precice.org/installation-overview.html)
- [rttr](https://www.rttr.org/)

## Building the adapter

For quick use, we provide a x86 shared library as download option. This `.so` file requires all the dependencies mentioned above. If you want to build the adapter yourself, checkout the FEBio repository and go to the branch of the version you want to build against. Set the cmake variable `FEBio_LIB_DIR` to this dir, as we need the include files to compile this plugin.
Then run:

```bash
cmake
make
```

Now there is a `libFEBioPreciceAdapter.so` in your cmake output folder. This is the file you need to refer to in your `febio.xml`.

## Using the adapter

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

### Configuration options

All attributes are contained in the `coupling_params` dict. The attributes are:

| Key | Possible Values | Description |
| --- | --- | --- |
| `element_set_to_couple` |  | Name of the element set that should be coupled via preCICE. Can be found in the febio input file. |
| `participant_name`  |  | Name of the preCICE participant of FEBio |
| `config_file_name`  |  |  Path to the preCICE configuration file  |
| `read_mesh_name` |  | Name of the mesh to read from |
| `read_data_name` |  | dict, see below |
| `write_mesh_name` |  | name of the mesh to write to |
| `write_data_name` |  | dict, see below |

#### write/read_data_name options

The adapter can read and write various data to and from FEBio. Possible data that can be handled are:

| Key | Possible Values | Description |
| --- | --- | --- |
| `type` | `function` `variable` | type of attribute for the FEMaterialPoint that is read or written. |
| `class_type` | `materialPoint` | This option sets what class type is accessed in FEBio, in future the possiblity of a `FEMaterial` or something like it could be added. |
| `febio_class_name` |  | Class name of the object that is being manipulated in FEBio e.g `FEFluidMaterialPoint` |
| `name` |  | Name of the attribute that is written/read `Vext` |
| `febio_type` | `vec3d` `vec2d` `double` `int` `float` `vector<dobule>` | data type of attribute being write/read |
| `mapping_name` |  | Name of data point in preCICE |
| `precice_type` | `scalar` `vector` | type of data point in precice |

### Example configuration

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

### Registering classes and generate.py

As we use reflection for the MaterialPoints, one has to register custom MaterialPoint classes with the plugin in `FEBioReflection.cpp` and extend two if statments. These statments are in `PreciceCallback::insertData` and `PreciceCallback::extractData`.

We added a `generate.py` that searches all subclasses of `FEMaterialPoint` in the current directory (your FEBio dir) and writes the needed code snippets to four files `FEBioReflection.cpp`, `includes` `insertData` and `extractData`. The `FEBioReflection.cpp` can directly be copied over the other three files need to be manually integrated into the `PreciceCallback.cpp`
To use `generate.py` you need to add the C++ header files to the include path

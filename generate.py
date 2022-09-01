#!/usr/bin/env nix-shell
#!nix-shell -i python -p python3Packages.pygccxml -p castxml -p llvm.dev -p libclang.dev -p libcxx.dev
from pygccxml import utils
from pygccxml import declarations
from pygccxml import parser
import os

# Find the location of the xml generator (castxml or gccxml)
generator_path, generator_name = utils.find_xml_generator()

# Configure the xml generator
xml_generator_config = parser.xml_generator_configuration_t(
            xml_generator_path=generator_path,cflags="-Wdelete-incomplete -Wpragma-once-outside-header -Wno-unused-command-line-argument",
                xml_generator=generator_name,
                include_paths=[Please add the path to your c++ headfiles])

def get_material_points_for_submodule(module_name):
    header_files = []
    header_files_path = []
    print(f"parsing {module_name}")
    for root, dirs, files in os.walk(module_name):
        for file in files:
            if file.endswith(".h") and file.startswith("FE") and not file.endswith("targetver.h") and not file.endswith("FEBioCommand.h"):
                filename = os.path.join(root,file)
                file_config = parser.file_configuration_t(
                    data=filename,
                    content_type=parser.CONTENT_TYPE.CACHED_SOURCE_FILE)
                header_files.append(file_config)
                header_files_path.append(filename)

    # Parse the c++ file
    decls = parser.parse(header_files, xml_generator_config,  compilation_mode=parser.COMPILATION_MODE.ALL_AT_ONCE)

    # Get access to the global namespace
    global_namespace = declarations.get_global_namespace(decls)
    material_points = {}
    for class_var in global_namespace.classes(allow_empty=True):
        for base in class_var.recursive_bases:
            if base.related_class.name == "FEMaterialPoint":
                material_points.update({class_var.name: class_var})
    return material_points, header_files_path

mps = {}
subfolders = [ f.path for f in os.scandir("./") if f.is_dir() ]
paths = []
for subfolder in subfolders:
    to_add, header_files_path = get_material_points_for_submodule(subfolder)
    paths = paths + header_files_path
    mps.update(to_add)
print(f"found {len(mps)} material point classes")
output_rttr_includes = '''
#include "FEBioReflection.h"
#include "PreciceCallback.h"
#include <rttr/registration>
'''

first = True
output_rttr = ''
output_insert = ''
output_extract = ''
output_includes = ''

output_rttr = output_rttr + 'RTTR_REGISTRATION {\n'
for key, value in mps.items():
    class_name = key
    output_includes = output_includes + f"#include <{value._location.file_name}>\n"
    if first:
        output_insert = output_insert + f'''
        if (className == "{class_name}") {{
            return wrapper.insert<{class_name}>(GetFEModel());
        }}
        '''
        output_extract = output_extract + f'''
        if (className == "{class_name}") {{
            return wrapper.extract<{class_name}>(GetFEModel());
        }}
        '''
        first = False
    else:
        output_insert = output_insert + f'''
        else if (className == "{class_name}") {{
            return wrapper.insert<{class_name}>(GetFEModel());
        }}
        '''
        output_extract = output_extract + f'''
        if (className == "{class_name}") {{
            return wrapper.extract<{class_name}>(GetFEModel());
        }}
        '''
    res = ''
    registration_start = f'registration::class_<{class_name}>("{class_name}")\n'
    res = res + registration_start
    variables = list(value.variables(allow_empty=True))
    for var_o in variables:
        var = var_o.name
        if var_o.access_type == "public":
            registration_prop = f'.property("{var}", &{class_name}::{var})\n'
            res = res + registration_prop

    functions = value.member_functions(allow_empty=True)
    for func_o in functions:
        method = func_o.name
        if func_o.access_type == "public":
            if len(func_o.overloads) > 0:
                print("We are currently not supporting overloaded functions")
            else:
                registration_method = f'.method("{method}", &{class_name}::{method})\n'
                res = res + registration_method
    res = res + ';\n'
    output_rttr = output_rttr + res
output_rttr = output_rttr + '}\n'
output_rttr = output_rttr_includes + output_includes + output_rttr
output_insert = output_insert + '''
    else {
        feLogError((std::string("Unknown Class Name ") + className).c_str());
    }
'''
text_file = open("FEBioReflection.cpp", "w")
n = text_file.write(output_rttr)
text_file.close()
text_file = open("insertData", "w")
n = text_file.write(output_insert)
text_file.close()
text_file = open("extractData", "w")
n = text_file.write(output_extract)
text_file.close()
text_file = open("includes", "w")
n = text_file.write(output_includes)
text_file.close()



'''
Author: git config user.name && git config user.email
Date: 2024-06-28 15:49:51
LastEditors: git config user.name && git config user.email
LastEditTime: 2024-09-09 15:42:33
FilePath: \tiger1.5\cmake\conanfile.py
Description: 

Copyright (c) 2024 by ${git_name_email}, All Rights Reserved. 
'''
from conan import ConanFile


class tigerRecipe(ConanFile):
    name = "tiger"
    description = "Required 3rd-libs for Zhe-Fong Grid Engine "
    settings = "os", "compiler", "build_type", "arch"
    generators = "CMakeToolchain", "CMakeDeps"
    build_policy = "missing"

    def requirements(self):
        # self.requires("eigen/3.4.0")
        # self.requires("vtk/9.2")
        # self.requires("spdlog/1.12.0")
        # self.requires("assimp/5.2.2")
        # self.requires("jsoncpp/1.9.5")
        # self.requires("glm/cci.20230113")
        # self.requires("gtest/1.14.0")
        # self.requires("opencascade/7.7.0")
        self.requires("cgns/3.3")
        self.requires("opencascade/7.7.0")
        # self.requires("autogrid/1.0")
        # self.requires("tigerlib/5.21-ubuntu")
        self.requires("libtiger/[>=9.1 <13.00]")

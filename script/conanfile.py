import os
from conan import ConanFile
from conan.tools.files import copy
from conan.tools.files import collect_libs
from datetime import datetime
class OpenCasCadeRecipe(ConanFile):
    name = "libtiger"
    version = ""
    settings = "os", "arch"
    build_policy = "never"
    package_type = "shared-library"

    def __init__(self,name):
        current_datetime = datetime.now()

        # 获取当前月份
        current_month = current_datetime.month

        # 获取当前日期
        current_day = current_datetime.day
        self.version=str(current_month)+"."+str(current_day)
        super().__init__()
        print(self.version)

        

    def layout(self):
        _os = str(self.settings.os).lower()
        _arch = str(self.settings.arch).lower()
        self.folders.build = os.path.join("libtiger", _os, _arch)
        self.folders.source = self.folders.build
        self.cpp.source.includedirs = ["include"]
        self.cpp.build.libdirs = ["lib"]
        self.cpp.build.bindirs = ["bin"]
        

    def package(self):
        local_include_folder = os.path.join(self.source_folder, self.cpp.source.includedirs[0])
        local_lib_folder = os.path.join(self.build_folder, self.cpp.build.libdirs[0])
        copy(self, "*", local_include_folder, os.path.join(self.package_folder, "include"), keep_path=True)
        copy(self, "*", local_lib_folder, os.path.join(self.package_folder, "lib"), keep_path=True)
        copy(self, "*", os.path.join(self.source_folder, "bin"), os.path.join(self.package_folder, "bin"), keep_path=True)
    def package_info(self):
        self.cpp_info.includedirs = ['include']  # Ordered list of include paths
        self.cpp_info.libdirs = ['lib']  # Directories where libraries can be found
        self.cpp_info.libs=collect_libs(self)
        self.cpp_info.bindirs = ['bin']  # Directories where executables and shared libs can be found
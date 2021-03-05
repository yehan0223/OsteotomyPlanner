cmake_minimum_required(VERSION 3.1)

project(OsteotomyPlanner)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/KitwareMedical/OsteotomyPlanner/wiki")
set(EXTENSION_CATEGORY "Osteotomy Planning")
set(EXTENSION_CONTRIBUTORS "Sam Horvath (Kitware Inc.), Johan Andruejol (Kitware Inc.)")
set(EXTENSION_DESCRIPTION "Osteotomy Planner is an extension designed to help the planning of bone surgeries.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/KitwareMedical/OsteotomyPlanner/master/Planner/Resources/Icons/Planner.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/KitwareMedical/OsteotomyPlanner/master/Screenshots/Screenshot.jpg")
set(EXTENSION_STATUS "In Development")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
# add_subdirectory(Planes)
# add_subdirectory(SplitModel)
# add_subdirectory(Planner)
add_subdirectory(ShrinkWrap)
add_subdirectory(OsteotomyModelToModelDistance)
add_subdirectory(ReplayPlan)
add_subdirectory(OsteotomyPlanner)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})

cmake_minimum_required(VERSION 3.10)
project(CylinderOrientation)

find_package(VTK REQUIRED)
vtk_module_config(VTK
  vtkCommonCore
  vtkIOImage
  vtkIOXML
  vtkInteractionWidgets
  vtkRenderingVolumeOpenGL2
)
include(${VTK_USE_FILE})

set(RENDERING_EXAMPLES_SRCS
  CylinderOrientation
)

foreach(name ${RENDERING_EXAMPLES_SRCS})
  add_executable(${name} MACOSX_BUNDLE ${name}.cxx)
  target_link_libraries(${name} ${VTK_LIBRARIES})
endforeach()
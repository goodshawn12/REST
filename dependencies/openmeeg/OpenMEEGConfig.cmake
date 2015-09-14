#-----------------------------------------------------------------------------
#
# OpenMEEGConfig.cmake -
# OpenMEEG CMake configuration file for external projects.
#
# This file is configured by OpenMEEG and used by the UseOpenMEEG.cmake module
# to load OpenMEEG's settings for an external project.

# The OpenMEEG source tree.
SET(OpenMEEG_SOURCE_DIR "/pipol/pipol/pipol7/openmeeg-release-2.1")

# The OpenMEEG include file directories.
SET(OpenMEEG_INCLUDE_DIRS "/usr/local/include/openmeeg")

# The OpenMEEG library directories.
# Includes OpenMEEG libs directory and OTHER libs directories
SET(OpenMEEG_LIBRARY_DIRS "/usr/local/lib")

# This includes all libraries that OpenMEEG generates.
SET(OpenMEEG_LIBRARIES "OpenMEEG;OpenMEEGMaths")

# The location of the UseOpenMEEG.cmake file.
SET(OpenMEEG_USE_FILE "/usr/local/lib/UseOpenMEEG.cmake")

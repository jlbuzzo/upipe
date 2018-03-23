# Main configuration file.

# Essential configuration paths.
SOURCES ?= ./src
MODULES ?= ./mdl
MODULE_NAME :=
LIBRARIES ?= ./lib
BINARIES ?= ./bin
DOCS ?= ./doc
LOGS ?= ./log


# Other conventions.
INPUT_DIR ?= ./inputs
OUTPUT_DIR ?= ./outputs


# Target chain variables.

# Target(s) for the verifycation and validation process of the inputs.
val_proc := val_proc

# Target(s) for the pre processes.
pre_proc := pre_proc

# Target(s) for the main processes.
mn_proc  := mn_proc

# Target(s) for the post processes.
pst_proc := pst_proc

# Target(s) for the extra processes.
ext_proc := ext_proc

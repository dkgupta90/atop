################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/TopologyOptimization/DensityValues.cpp \
../source/TopologyOptimization/RefinementIndicator.cpp \
../source/TopologyOptimization/cell_prop.cpp \
../source/TopologyOptimization/design_analysis.cpp \
../source/TopologyOptimization/neighbors.cpp \
../source/TopologyOptimization/penalization.cpp 

OBJS += \
./source/TopologyOptimization/DensityValues.o \
./source/TopologyOptimization/RefinementIndicator.o \
./source/TopologyOptimization/cell_prop.o \
./source/TopologyOptimization/design_analysis.o \
./source/TopologyOptimization/neighbors.o \
./source/TopologyOptimization/penalization.o 

CPP_DEPS += \
./source/TopologyOptimization/DensityValues.d \
./source/TopologyOptimization/RefinementIndicator.d \
./source/TopologyOptimization/cell_prop.d \
./source/TopologyOptimization/design_analysis.d \
./source/TopologyOptimization/neighbors.d \
./source/TopologyOptimization/penalization.d 


# Each subdirectory must supply rules for building sources it contributes
source/TopologyOptimization/%.o: ../source/TopologyOptimization/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/dkgupta/bin/deal.II/include -I"/home/dkgupta/WORK/projects/atop/backupproj/include" -O0 -g3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



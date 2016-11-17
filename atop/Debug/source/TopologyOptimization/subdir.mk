################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/TopologyOptimization/DensityValues.cpp \
../source/TopologyOptimization/RefinementIndicator.cpp \
../source/TopologyOptimization/adaptivity.cpp \
../source/TopologyOptimization/cell_prop.cpp \
../source/TopologyOptimization/designField.cpp \
../source/TopologyOptimization/design_analysis.cpp \
../source/TopologyOptimization/neighbors.cpp \
../source/TopologyOptimization/optimizedesign.cpp \
../source/TopologyOptimization/penalization.cpp \
../source/TopologyOptimization/projection.cpp 

OBJS += \
./source/TopologyOptimization/DensityValues.o \
./source/TopologyOptimization/RefinementIndicator.o \
./source/TopologyOptimization/adaptivity.o \
./source/TopologyOptimization/cell_prop.o \
./source/TopologyOptimization/designField.o \
./source/TopologyOptimization/design_analysis.o \
./source/TopologyOptimization/neighbors.o \
./source/TopologyOptimization/optimizedesign.o \
./source/TopologyOptimization/penalization.o \
./source/TopologyOptimization/projection.o 

CPP_DEPS += \
./source/TopologyOptimization/DensityValues.d \
./source/TopologyOptimization/RefinementIndicator.d \
./source/TopologyOptimization/adaptivity.d \
./source/TopologyOptimization/cell_prop.d \
./source/TopologyOptimization/designField.d \
./source/TopologyOptimization/design_analysis.d \
./source/TopologyOptimization/neighbors.d \
./source/TopologyOptimization/optimizedesign.d \
./source/TopologyOptimization/penalization.d \
./source/TopologyOptimization/projection.d 


# Each subdirectory must supply rules for building sources it contributes
source/TopologyOptimization/%.o: ../source/TopologyOptimization/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++-5 -DDEAL_II_USE_CXX11 -I/home/dkgupta/WORK/projects/atop/atop/include -I/home/dkgupta/bin/deal.II/include -I/usr/include/suitesparse -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



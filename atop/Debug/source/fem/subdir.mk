################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/fem/boundary_values.cpp \
../source/fem/create_design.cpp \
../source/fem/define_mesh.cpp \
../source/fem/fem.cpp \
../source/fem/output.cpp 

OBJS += \
./source/fem/boundary_values.o \
./source/fem/create_design.o \
./source/fem/define_mesh.o \
./source/fem/fem.o \
./source/fem/output.o 

CPP_DEPS += \
./source/fem/boundary_values.d \
./source/fem/create_design.d \
./source/fem/define_mesh.d \
./source/fem/fem.d \
./source/fem/output.d 


# Each subdirectory must supply rules for building sources it contributes
source/fem/%.o: ../source/fem/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++-5 -DDEAL_II_USE_CXX11 -I/home/dkgupta/WORK/projects/atop/atop/include -I/home/dkgupta/bin/deal.II/include -I/usr/include/suitesparse -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



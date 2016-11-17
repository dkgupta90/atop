################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/examples/step1.cpp 

OBJS += \
./source/examples/step1.o 

CPP_DEPS += \
./source/examples/step1.d 


# Each subdirectory must supply rules for building sources it contributes
source/examples/%.o: ../source/examples/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/dkgupta/bin/deal.II/include -I"/home/dkgupta/WORK/projects/atop/backupproj/include" -O0 -g3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



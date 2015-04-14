################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../examples/simpleTO/step1.cpp 

OBJS += \
./examples/simpleTO/step1.o 

CPP_DEPS += \
./examples/simpleTO/step1.d 


# Each subdirectory must supply rules for building sources it contributes
examples/simpleTO/%.o: ../examples/simpleTO/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/dkgupta/bin/deal.II/include -I/home/dkgupta/WORK/projects/mtop/topopt/include -O0 -g3 -Wall -c -fmessage-length=0 -std=gnu++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



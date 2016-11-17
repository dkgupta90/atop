################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/data_format/bitmap.cpp 

OBJS += \
./source/data_format/bitmap.o 

CPP_DEPS += \
./source/data_format/bitmap.d 


# Each subdirectory must supply rules for building sources it contributes
source/data_format/%.o: ../source/data_format/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/home/dkgupta/bin/deal.II/include -I"/home/dkgupta/WORK/projects/atop/backupproj/include" -O0 -g3 -Wall -c -fmessage-length=0 -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/fem/ErrorIndicator/stressJumpIndicator.cpp 

OBJS += \
./source/fem/ErrorIndicator/stressJumpIndicator.o 

CPP_DEPS += \
./source/fem/ErrorIndicator/stressJumpIndicator.d 


# Each subdirectory must supply rules for building sources it contributes
source/fem/ErrorIndicator/%.o: ../source/fem/ErrorIndicator/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++-5 -DDEAL_II_USE_CXX11 -I/home/dkgupta/WORK/projects/atop/atop/include -I/home/dkgupta/bin/deal.II/include -I/usr/include/suitesparse -I/home/dkgupta/IpOpt/CoinIpopt/build/include -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



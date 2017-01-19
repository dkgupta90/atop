################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/optimizer/ipopt_interface.cpp \
../source/optimizer/optimality_criteria.cpp 

OBJS += \
./source/optimizer/ipopt_interface.o \
./source/optimizer/optimality_criteria.o 

CPP_DEPS += \
./source/optimizer/ipopt_interface.d \
./source/optimizer/optimality_criteria.d 


# Each subdirectory must supply rules for building sources it contributes
source/optimizer/%.o: ../source/optimizer/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++-5 -DDEAL_II_USE_CXX11 -I/home/dkgupta/WORK/projects/atop/atop/include -I/home/dkgupta/bin/deal.II/include -I/usr/include/suitesparse -I/home/dkgupta/IpOpt/CoinIpopt/build/include -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



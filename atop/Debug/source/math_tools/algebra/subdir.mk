################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/math_tools/algebra/MatrixVector.cpp \
../source/math_tools/algebra/integration.cpp 

OBJS += \
./source/math_tools/algebra/MatrixVector.o \
./source/math_tools/algebra/integration.o 

CPP_DEPS += \
./source/math_tools/algebra/MatrixVector.d \
./source/math_tools/algebra/integration.d 


# Each subdirectory must supply rules for building sources it contributes
source/math_tools/algebra/%.o: ../source/math_tools/algebra/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++-5 -DDEAL_II_USE_CXX11 -I/home/dkgupta/WORK/projects/atop/atop/include -I/home/dkgupta/bin/deal.II/include -I/usr/include/suitesparse -I/home/dkgupta/IpOpt/CoinIpopt/build/include -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



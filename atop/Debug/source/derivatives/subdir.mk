################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../source/derivatives/compliance.cpp \
../source/derivatives/compliant_mechanism.cpp \
../source/derivatives/sensitivity.cpp 

OBJS += \
./source/derivatives/compliance.o \
./source/derivatives/compliant_mechanism.o \
./source/derivatives/sensitivity.o 

CPP_DEPS += \
./source/derivatives/compliance.d \
./source/derivatives/compliant_mechanism.d \
./source/derivatives/sensitivity.d 


# Each subdirectory must supply rules for building sources it contributes
source/derivatives/%.o: ../source/derivatives/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++-5 -DDEAL_II_USE_CXX11 -I/home/dkgupta/WORK/projects/atop/atop/include -I/home/dkgupta/bin/deal.II/include -I/usr/include/suitesparse -I/home/dkgupta/IpOpt/CoinIpopt/build/include -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '




#This is for all systems
limit stacksize unlimited

# This is for IBM AIX systems
setenv XLSMPOPTS "stack=20000000"

# This is for Linux systems 
setenv KMP_STACKSIZE 20000000

# This is for HP/Compaq Tru64 systems
setenv MP_STACK_SIZE 20000000

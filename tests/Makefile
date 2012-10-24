all: errmod_test.o
	gcc -Wall -pedantic -o errmod_test errmod_test.o ../errmod.o -lm
errmod_test.o:
	gcc -c errmod_test.c -I ..

realmain:main_iplib.c ip_lib.o bmp.o
	gcc main_iplib.c ip_lib.o bmp.o -o realmain --ansi -pedantic -ggdb -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra 

ip_lib.o: ip_lib.c ip_lib.h bmp.o
	gcc ip_lib.c -o ip_lib.o --ansi -pedantic -Wall -c -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra 
	
bmp.o: bmp.c bmp.h
	gcc bmp.c -o bmp.o -Wall -c -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra 

clean: 
	@rm -f bmp.o mainbmp ip_lib.o main primaprovameta.bmp

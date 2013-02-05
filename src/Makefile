PREFIX			= /usr/local/bin
TEST				= test/
CC					= gcc  
CFLAGS			= -O9 -Wall -m64 -std=gnu99
LD					= -lpthread -lm
DEBUGFLAGS		= -ggdb -Wall -m64 -std=gnu99

SOURCES			= main.c interface.c nuc_code.c misc.c io.c

PROGS				= tagdust
DEBUGPROGS	= tagdust_debug
OBJECTS			= $(SOURCES:.c=.o)
DEBUGOBJECTS	= $(SOURCES:.c=_debug.o)

MANDIR      = /usr/share/man/
MANSUFFIX = 1

all: $(OBJECTS)
	$(CC) $(CFLAGS) $(LD) $(OBJECTS) -o $(PROGS)
	
%.o:	%.c
	$(CC) $(CFLAGS) -c $< -o $@

debug: $(DEBUGOBJECTS)
	$(CC) $(DEBUGFLAGS) $(LD) $(DEBUGOBJECTS) -o $(DEBUGPROGS)	

%_debug.o:	%.c
	$(CC) $(DEBUGFLAGS) -c $< -o $@
	

install:
	cp delve.1 $(MANDIR)/man$(MANSUFFIX)/delve.1
	gzip $(MANDIR)/man1/delve.1
	cp $(PROGS) $(PREFIX)
	

clean: 
	@\rm -f $(PROGS) $(OBJECTS)
	@\rm -f $(DEBUGPROGS) $(DEBUGOBJECTS)
	@\rm -f $(THREADPROGS) $(THREADOBJECTS)
	@\rm -f *~
	

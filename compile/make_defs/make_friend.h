#=================================================================
#               FRIEND_FILES
#=================================================================


#=================================================================
cmalloc.o     :          $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(CODE)/friend_lib/cmalloc.c
	$(ECHO) $@
	$(COBJ) $(CODE)/friend_lib/cmalloc.c

#-----------------------------------------------------------------
friend_lib.o   :         $(STANDARD) $(DEFINES) \
                         $(FRND_ENT) \
                         $(CODE)/friend_lib/friend_lib.c
	$(ECHO) $@
	$(COBJ) $(CODE)/friend_lib/friend_lib.c

#=================================================================

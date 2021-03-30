#!/usr/bin/env bash

#this generates dependency files for a fortran file to stdout
# $1 = source file
# $2 = root build directory (can be blank, a value or a make variable unresolved)
# $3 - if present only use modules (*.mod) files as dependencies, otherwise use both *.mod and its corresponding object (*.o) file as dependencies

#for the given source file, 
#  1) for each module "use"d, if source file (*.f) that defines the module exists then
#       a) if $3 is present, add the corresponding object file that defines the module to dependency list
#              this will force more recompiling (just basing off object or *.mod files is not enough for general fortran compilers)
#              if using gfortran >= 4.3.0 then *.mod timestamp is only changed if the module interface changed, 
#                therefore we can avoid including the *.o file(s) in dependencies and simply use the *.mod files.
#                in this case you should pass in anything for $3
#                note, this is NOT the behavior of intel's ifort
#       b) add the *.mod file to dependency list
#  2) for each fortran 'include', add that also to the dependencies list (if it exists) and recursively do the dependencies of that included file
#  3) for each module defined in the source file, add a simple rule for *.mod with dependency of *.o file

# NOTE: 
#  we don't deal with preprocessing or preprocessing includes (i.e. "#include ...')
#  all *.mod files are placed in $2 (= "$(BLDDIR)/", say) and all *.o are mirrored w.r.t. *.f (e.g. src/plot/poly.f and $(BLDDIR)/src/plot/poly.o) 

find_mod_inc() {

	deps=""

	#loop over each module 'use'd in this source file: $j will be the name (lower case) of the module
	for j in $(grep -i "^ *use "  $1 | awk '{print tolower($0)}' | sed "s/,/ /g; s/::/ /; s/ non_intrinsic / /; s/ intrinsic / /" | \
		awk '{print $2}' | sort -u | sed "s/\!/ /" | awk '{print $1}'|sort -u); do

        	#for each module found, find the source file that defines it- may not exist (or may be in the current source file)
		moddef=""
       	 	for i in $(find . -name "*.f" | grep -v $1 | sed "s/^\.\///" | xargs grep -i "^ *module *${j}" | \
			awk -F: '{print $1}'|sort -u); do   #i = name of source file w/ possible module definition
			(( $(grep -i "^ *module *${j}" $i | awk '{print tolower($2)}' | \
				sort -u | sed "s/\!/ /" | awk -v m=$j '{if ($1 == m) print m}' | wc -l) > 0 )) && \
				moddef="1" && [ ! -n "${onlymod}" ] && deps+=" ${blddir}${i/.f/.o}"  #add correspond source object to dependencies
        	done

        	#add *.mod file to dependency if module is not defined in the source file
		[ -n "$moddef" ] && (( $(grep -i "^ *module *${j}" $1 | awk '{print tolower($2)}' | \
			sort -u | sed "s/\!/ /" | awk -v m=$j '{if ($1 == m) print m}' |wc -l) < 1 )) && \
			deps+=" ${blddir}${j}.mod"

	done

	#look for any 'include' statements... and add them onto the dependencies if they are in any subdirectories
	for i in $(grep -i "^ *include " $1 | sed "s/\"/ /g; s/\'/ /g"| awk '{print $2}'| sort -u);do
		j=$(basename $i) && (( $(find . -name $j | wc -l) > 0 )) && deps+=" $(dirname $1)/${i}$(find_mod_inc $(dirname $1)/${i})" #recursion!
	done
	echo "${deps}"
}

[ -n "$3" ] && onlymod="1"
blddir="$2" && sfile="$1" && mkline="${blddir}${sfile/.f/.o} : ${sfile}$(find_mod_inc ${sfile})" #object file : source file
echo "$mkline"  #ok, final dependency list for the object file

#for each module defined in the source file, add a simple rule for *.mod w/ dependency of the object file
for i in $(grep -i "^ *module " $sfile | awk '{if (tolower($1) == "module" && tolower($2) != "procedure") print tolower($2)}' | \
	sort -u | sed "s/\!/ /" | awk '{print $1}');do 
	echo && echo "${blddir}${i}.mod : ${blddir}${sfile/.f/.o}"
done

Κ23γ: Ανάπτυξη Λογισμικού για Αλγοριθμικά Προβλήματα Χειμερινό εξάμηνο 2022-23 


2η Προγραμματιστική Εργασία Πολυγωνοποίηση σημειοσυνόλου με τη χρήση της βιβλιοθήκης CGAL (C++)


Παναγιώτης Κοντοειδής 1115201900266
Στέλιος Δημητριάδης 1115201900050


Σκοπός της εργασίας είναι η πολυγωνοποίηση βέλτιστης επιφάνειας σημειοσυνόλου S το οποίο περιέχει n σημεία στο χώρο.


Τα αρχεία που περιλαμβάνονται είναι τα: 


PROJECT_2.cpp
CMakeLists.txt
output.txt
output_polygon.txt
output_polygon_subdiv.txt
uniform-0000xx-1.instance
include/Area_maximization_minimization.cpp
include/Area_maximization_minimization.hpp
include/CMakeLists.txt
README.txt
README.pdf
.git




















Οδηγίες μεταγλώττισης του προγράμματος:


1. 
cmake -DCGAL_DIR=$CMAKE_INSTALLED_PREFIX/lib/CGAL -DCMAKE_BUILD_TYPE=Release .


2.        
make


3.
 ./PROJECT_2 -i uniform-0000100-1.instance(input must have this kind of name ex xxx-0000zzz-x.instance zzz=number of points) -o output.txt -algorithm (simulated_annealing or local_search) -L 500(number defined according to algorithm) -max(-max or -min) ((-threshold <double>(for local search)) or -annealing <>(local global or subdivision)) ./PROJECT -i uniform-0000100-1.instance(same input as PROJECT_2) -o output_polygon.txt(this must be strictly this name) -algorithm <>(incremental or convex_hull) -edge_selection <>(1 or 2 or 3) -init <>(1a or 2a or 1b or 2b)




Στους παρακάτω πίνακες αναφέρονται οι μέσοι χρόνοι εκτέλεσης και το ratio των δοθείσων παραδειγμάτων. Τα αποτελέσματα προέκυψαν από την επαναλαμβανόμενη εκτέλεση των αλγορίθμων. Τα συμπεράσματα που προκύπτουν είναι τα εξης:




Για τον αλγόριθμο Προσομοιωμένης Ανόπτησης με annealing local, με την δημιουργια πολυγωνου με αλγοριθμο convex hull, με random επιλογη σημειων, στο μικρο file το ratio μικραίνει κατα πολυ, με min το ratio μεγαλώνει μαζί και ο χρονος ενω, με max εχουμε παρομοια συμπεριφορα με το random. Με την χρηση incremental οι χρονοι ειναι παρόμιοι ενώ το ratio με random δεν διαφέρει απο το προηγουμενο οπως με convex hull, τα υπολοιπα παραμενουν παρόμοια. Για max με convex hull δεν υπάρχει μεγάλη διαφορα ratio όσο με min και παρατηρω παλι οτι με max εχω μικροτερο ratio απο το αρχικο. Με annealing global ο χρόνος μεγαλώνει αρκετά και στα δυο files ομως τα αποτελεσματα ειναι αρκετα καλυτερα σε ολες τις περιπτώσεις ασχέτως τον αλγόριθμο convex hull ή incremental και ασχετως την min max random επιλογή. Σχετικά με το μεγαλο αρχειο με global δεν υπάρχουν πολλα δεδομενα καθως κανει αρκετη ωρα να τρεξει. Συμπεραίνουμε πως με local annealing υπάρχει πολύ μεγαλύτερη αποδοτικότητα απο οτι με global και οι διαφορες του local με global ratio ειναι μικρες αλλα καλύτερες στο global.




Για τον αλγόριθμο  Τοπικής Αναζήτησης, με την δημιουργία πολυγώνου με τον αλγόριθμο convex hull, στο μικρό file ratio παρατηρούμε ότι ο χρόνος σε σύγκριση με τον simulated annealing local είναι ελάχιστος ενώ σε αρχεία με μεγάλο file ratio ο χρόνος αυξάνεται κατά πολύ, ξεπερνώντας τον όχι όμως όταν είναι global. Αυτό μπορεί να δικαιολογηθεί στην πολυπλοκότητα της Τοπικής αναζήτησης όπου στην υλοποίηση μας είναι Ο(ν^4).Οι επιλογές ταξινόμησης των σημείων δεν προσφέρουν κάποια αξιοσημείωτη διαφορά μεταξύ τους στα τελικά αποτελέσματα για αυτό για την δειγματοληψία παρουσιάζεται 1α επιλογή της πρώτης εργασίας, δηλαδή η φθίνουσα ταξινόμηση κατά χ.








ΑΛΓΌΡΙΘΜΟΣ:SIMULATED-ANNEALING
LOCAL
L=5000
FILE:
	euro-night-0000030.instance
	uniform-0000900-2.instance
	MIN:
CONVEX_HULL:
RANDOM


	87 ms
ratio initial: 0.615098
ratio new: 0.35112
	4343 ms
ratio initial: 0.509701
ratio new: 0.442568
	MIN:
CONVEX_HULL:
MIN


	173 ms
ratio initial: 0.590296
ratio new: 0.644699
	4904 ms
ratio initial: 0.391436
ratio new: 0.384805
	MIN:
CONVEX_HULL:
MAX
	115 ms
ratio initial: 0.783189
ratio new: 0.407124
	4188 ms
ratio initial: 0.6911
ratio new: 0.519274
	MIN:
INCREMENTAL:
RANDOM
INIT:1a
	95 ms
ratio initial: 0.537495
ratio new: 0.400635
	4371 ms
ratio initial: 0.508887
ratio new: 0.428958
	MIN:
INCREMENTAL:
MIN
INIT:2a
	106 ms
ratio initial: 0.27611
ratio new: 0.285075
	5006 ms
ratio initial: 0.391436
ratio new: 0.379648
	MIN:
INCREMENTAL:
MAX
INIT:2b
	120 ms
ratio initial: 0.802732
ratio new: 0.516447
	4049 ms
ratio initial: 0.6911
ratio new: 0.508905
	MAX:
CONVEX_HULL:
RANDOM
	104 ms
ratio initial: 0.636084
ratio new: 0.701211
	3654 ms
ratio initial: 0.513851
ratio new: 0.56701
	MAX:
CONVEX_HULL:
MIN
	159 ms
ratio initial: 0.590296
ratio new: 0.704599
	4472 ms
ratio initial: 0.391436
ratio new: 0.510648
	MAX:
CONVEX_HULL:
MAX
	131 ms
ratio initial: 0.783189
ratio new: 0.705141
	3768 ms
ratio initial: 0.6911
ratio new: 0.640461
	MAX:
INCREMENTAL:
RANDOM
INIT:1a
	99 ms
ratio initial: 0.43299
ratio new: 0.584701
	3513 ms
ratio initial: 0.503059
ratio new: 0.567066
	MAX:
INCREMENTAL:
MIN
INIT:2a
	103 ms
ratio initial: 0.27611
ratio new: 0.470112
	3428 ms
ratio initial: 0.62563
ratio: 0.6301
	MAX:
INCREMENTAL:
MAX
INIT:2b
	168 ms
ratio initial: 0.802732
ratio new: 0.74848
	3696 ms
ratio initial: 0.615842
ratio new: 0.611745
	



ΑΛΓΌΡΙΘΜΟΣ:SIMULATED-ANNEALING
GLOBAL
L=5000


FILE:
	euro-night-0000030.instance
	uniform-0000900-2.instance
	MIN:
CONVEX_HULL:
RANDOM
	560 ms
ratio initial: 0.596759
ratio new: 0.326261
	297296 ms
ratio initial: 0.502292
ratio new: 0.347179
	MIN:
CONVEX_HULL:
MIN
	537 ms
ratio initial: 0.590296
ratio new: 0.383739
	

	MIN:
CONVEX_HULL:
MAX
	572 ms
ratio initial: 0.783189
ratio new: 0.388654
	

	MIN:
INCREMENTAL:
RANDOM
INIT:1a
	508 ms
ratio initial: 0.636522
ratio new: 0.447771
	

	MIN:
INCREMENTAL:
MIN
INIT:2a
	511 ms
ratio initial: 0.802732
ratio: 0.430266
	

	MIN:
INCREMENTAL:
MAX
INIT:2B
	578 ms
ratio initial: 0.27611
ratio new: 0.317933
	

	MAX:
CONVEX_HULL:
RANDOM
	510 ms
ratio initial: 0.619476
ratio new: 0.672609
	375346 ms
ratio initial: 0.500049
ratio new: 0.913725
	MAX:
CONVEX_HULL:
MIN
	517 ms
ratio initial: 0.590296
ratio new: 0.726048
	

	MAX:
CONVEX_HULL:
MAX
	593 ms
ratio initial: 0.783189
ratio new: 0.671705
	

	MAX:
INCREMENTAL:
RANDOM
INIT:1a
	602 ms
ratio initial: 0.658003
ratio new: 0.712046
	

	MAX:
INCREMENTAL:
MIN
INIT:2a
	571 ms
ratio initial: 0.802732
ratio new: 0.6359
	

	MAX:
INCREMENTAL:
MAX
INIT:2b
	596 ms
ratio initial: 0.27611
ratio new: 0.77166
	



	



ΑΛΓΟΡΙΘΜΟΣ: LOCAL SEARCH


L=1
THRESHOLD =0.1
MAX


FILE:
	uniform-0000015-1.instance
	uniform-0000900-2.instance
	MIN:CONVEX_HULL:RANDOM
	1ms
ratio initial:0.738381
ratio new:0.715806
	21174ms
ratio initial:0.494375
ratio new:0.494174
	MIN:CONVEX_HULL:MIN
	1ms
ratio initial:0.738381
ratio new:0.715806
	19956ms
ratio initial:0.391436
ratio new:0.391267
	MIN:CONVEX_HULL:
MAX
	1ms
ratio initial:0.738381
ratio new:0.715806
	19846ms
ratio initial:0.6911
ratio new:0.70245
	MIN:INCREMENTAL:RANDOM 1a
	1ms
ratio initial:0.738381
ratio new:0.715806
	19033ms
ratio initial:0.496506
ratio new:0.49684
	MIN:INCREMENTAL:MIN 1a
	1ms
ratio initial:0.738381
ratio new:0.715806
	19224ms
ratio initial:0.60994
ratio new:0.60885
	MIN:INCREMENTAL:MAX 1a
	1ms
ratio initial:0.738381
ratio new:0.715806
	23556ms
ratio initial:0.347161
ratio new:0.347494

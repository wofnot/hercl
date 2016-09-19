------ README ------

---------- Introduction ----------

HERCL is a programming language designed for evolutionary computation with the specific aim of allowing new programs to be created by combining patches of code from different parts of other programs, at multiple scales.

HERCsearch is a global random search strategy, based on Linear Genetic Programming, used to evolve HERCL programs to perform specific tasks.

---------- Supervised Learning Mode ----------

In supervised learning, agents (HERCL programs) are given training and test data. Agents are evaluated and evolved based on their ability to produce the correct output when given the corresponding input, as these are defined by the training data. The test data does not affect evolution but is used to judge how succesfully the process produced a general solution to the task.

---------- Interactive Learning Mode ----------

In interactive learning, agents interact with an evironment (usually written in C or C++), which evaluates them based on criteria defined in the environment. This is to allow evolution on tasks that can't be represented easily as input/output pairs (i.e. what supervised learning does). Examples of tasks that are suited to interactive environments are the double-pole balancing problem, or exploring a maze.

---------- Additional Information ----------

See the other .txt files in /doc for further information on HERCsearch and HERCL.

http://www.cse.unsw.edu.au/~blair/pubs/2013BlairCEC.pdf
http://www.cse.unsw.edu.au/~blair/pubs/2014BlairGECCO.pdf
http://www.cse.unsw.edu.au/~blair/pubs/2015BlairACALCI.pdf

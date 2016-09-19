----- About -----

(Skip to 'How to' below if all you want to do is run a co-evolution sequence)

Alternately evolves artists and art critics (both are HERCL programs). The critics are trained to distinguish 'real' images from HERCL-generated ones, while the artists are evaluated based on their ability to 'fool' critics. Running this process iteratively creates a competitive dynamic between the artist and the critic, forcing the artist to generate novel images at each phase. Based on an approach described by Machado & Romero in "On Evolutionary Computer-Generated Art" (2011).

The details of the process are as follows:

1. We start with a set of positive items and a set of negative items. Positive items are images we would like to be considered 'real' images by the critic. Negative items are the opposite. In our experiments so far we use 100-200 images taken from the internet as positive examples, and a blank image as a single negative example.

2. Each image is condensed into a feature vector based on various image metrics (see src/metrics/metrics.h for more info). This is process is done by make-train.sh.

3. With the positive and negative feature vectors as training data, supervised learning is used to evolve a HERCL program (the 'critic') to classify images as either 'real' or 'fake'. Usually this evolution is run several times with different a different random seed, giving a diverse set of critics. make-critics.sh handles this step.

4. Now we evolve another HERCL program (the 'artist') whose goal is to produce images that will fool the critics that have just been evolved. When an artist is being evaluated, its output (in the form of simple numerical messages) is interpreted as a sequence of instructions in a simple 'turtle' graphics environment (see src/inter/art-school.h for more info). When the artist has finished, the produced image is converted to a feature vector and given to the previously evolved critics. The average score given by the critic determines the fitness function for the artist. Once the artist is able to produce an image that is rated highly enough by critics, it is considered successful and the image is rendered and saved. Again, we run the evolution multiple times with different seeds to produce a set of artists. make-artists.sh does this bit.

(Steps 2-4 describe a single 'phase'. Usually the results after 1 phase will be uninteresting - it is not difficult for the critic to distinguish a blank image from a real image, and so the critic's logic can be very simple, making it easy for the artist to fool. Consequently the drawings produced in phase 1 are often very simple (single dots or lines).)

5. Now, after the completion of the phase, we take the generated art and add it to the _negative_ dataset for the next phase, and repeat steps 2-4. This means the new critic must be more sophisticated in order to distinguish the real images from fake ones. Over several phases, both the artist's and the critic's task become more difficult, resulting in increasingly complex images.

----- How to -----

To run a coevolution sequence:

0. Compile the necessary executables: runMetrics, hercsearch and art-school.

1. First you need to create some initial training data
  (NB Currently all images should be 128x128)
  a) place positive examples in data/<dataset name>/pos/
  b) place negative examples in data/<dataset name>/neg/ (just a single white image works fine)
  
2. Create a folder for the sequence to be generated in.
  a) Within this folder create a file settings.txt (see settings.template for structure)
     Importantly, the field train_data should refer to the folder created in step 1 (data/<dataset name>)
3. Run a phase
  a) Run ./run-phase.sh <folder> 1
  b) (to run a number of phases in sequence)
   for i in {1..20}; do ./run-phase.sh folder $i; done;

 (NB If a phase is interrupted midway, it should be fine to resume it, but you must find the latest file,
  e.g <folder>/critics/5.hrc or <folder>/artists/3.hrc, which should be empty, and remove it.)

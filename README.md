# Weighted Correlation Atom Decomposition 
October 2015

WCAD is an intonation modeling algorithm based on using the weighted correlation as a cost function in the matching pursuit framework. 

To run WCAD, clone it from `git` and start `wcad.py`.

```sh
git clone git@github.com:dipteam/wcad.git
cd wcad
# Edit paths.cfg to match your configuration
./wcad.py $path_to_audio/$my_wav.wav $path_to_results
```
To run WCAD you need a compiled binary for the pitch tracker `compute-kaldi-pitch-feats` from [Kaldi](http://kaldi-asr.org/). WCAD uses it to get a continuous pitch estimate and a probability of voicing. Kaldi's code can be found [here](https://github.com/kaldi-asr/kaldi). Once compiled you need to add the path to the binary in the file `./paths.cfg`.

The code supports [Praat](http://www.fon.hum.uva.nl/praat/) TextGrid annotations for extracting the utterance boundaries and for plotting syllables/words/phones in the final plots. By default these should be in a folder called `annotations/` in the audio folder. This functionality is added with the inclusion of `textgrid.py` from the [Natural Language Toolkit (NLTK)](https://github.com/nltk/nltk_contrib), which can be found [here](https://github.com/nltk/nltk_contrib/blob/master/nltk_contrib/textgrid.py).  By default, the use of TextGrid files is deactivated.  To use TextGrid functionalities, please edit `wcad/object_types/params.py` by setting `self.plot_textgrids` to True.  You can also set the speech start and end detection used for phrase extraction to use the TextGrid annotations instead of energy, by setting `self.fix_pos_ref_func` to `'energy'`.


Most of the algorithm parameters can be tweaked in the file `wcad/object_types/params.py`.

The Weighted Correlation Atom Decomposition algorithm is described in the paper:

Branislav Gerazov, Pierre-Edouard Honnet, Aleksandar Gjoreski, and Philip N. 
Garner, "Weighted Correlation based Atom Decomposition Intonation Modelling," in Proc. of Interspeech 2015, Dresden, Germany, Sep 06 - 10, 2015.

```
@inproceedings{Gerazov2015,
    author = {Gerazov, Branislav and Honnet, Pierre-Edouard and Gjoreski, Aleksandar and Garner, Philip N.},
    title = {Weighted Correlation based Atom Decomposition Intonation Modelling},
    booktitle = {Proceedings of Interspeech},
    location = {Dresden, Germany},
    month = {Sep},
    year = {2015}
}
```

Branislav Gerazov

[DIPteam](http://dipteam.feit.ukim.edu.mk/)

[Faculty of Electrical Engineering and Information Technologies](http://feit.ukim.edu.mk)

Ss Cyril and Methodius University of Skopje,

Macedonia

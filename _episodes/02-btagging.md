---
title: "B-Tagging Exercise"
teaching: 15 min
exercises: 10 min
questions: "How are b hadrons discriminated in analyses"
-
.objectives:
-"Understand what b tagging is"
-"Learn b tagging algorithms"
-"be able to access and make cuts on b-tag discriminators"

.keypoints:
-"b tagging algorithms discriminate b jets from jets produced by the hadronization of light quarks and gluon"
-"b tagging algorithms produce a disriminator value for each jet which measures the likelihood of a b jet"
-"The BTV POG defines operating working points for the b-tagging algorithms to define cuts in the discriminator and
are based on a misidentification probability for light-flavor jets"
---



##
Jet reconstruction and identification is an important part of the analyses at the LHC.
Several b tag algorithms exist all with the goal of idetifying jets from the hadronization of b-quarks. 
B jets must be discriminated from jets produced by the hadronization of light quarks and gluon.
The algorithms first associate good quality tracks (vector of references to Tracks) with jets.
Tracks are associated with jets using a cone-based association or the explicit jet-track association. 
This is called jet-track association (JTA). Track quality criteria are then applied and tracks are passed through 
track-based tagging algorithms. If secondary vertices are reconstructed, they can be passed to SV-based tagging algorithms or to combined tagging algorithms. 
They produce a single, real number called a btag "discriminator" for each jet. 
If the discriminator is more positive then the higher the b-jet probability. 

# B Tagging Algorithms

The specific details depend upon the algorithm use. However, they all exploit a property (or properties) of b hadrons namely:

* the long lifetime
* large mass
* high track multiplicity
* large semileptonic branching fraction
* hard fragmentation fuction 

Algorithms that are used for b-tagging:

*Track Counting - This is a simple tag that identifies a b jet if it contains at least N tracks with significance of the impact parameter exceeding some value

*Jet Probability - This combines information from all selected tracks in the jet and uses probability density functions to assign a probability to each track

*Soft Muon and Soft Electron- This identifies a b by searching for a lepton from the semi-leptonic b decay mode

*Simple Secondary Vertex- This reconstructs the b decay vertex and uses kinematic variables related to it to 
calculate a discriminator

*Combined Secondary Vertex- This is a complex tag that exploits all known kinematic variables of the b jets, information
about impact parameter significance and the secondary vertex to distinguish b jets.

# Accessing b-tag information from NanoAOD

The following lines are taken from the [AOD2NanoAODOutreachTool](https://github.com/cms-opendata-analyses/AOD2NanoAODOutreachTool) and show how to access b-tag information.

~~~
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"

//define jet handle
Handle<CaloJetCollection> jets;

//get jets from the event
iEvent.getByLabel(InputTag("ak5CaloJets"), jets);

//define b-tag discriminators handle
Handle<JetTagCollection> btags;

//get b-tag discriminators from the event
iEvent.getByLabel(InputTag("combinedSecondaryVertexBJetTags"), btags);

//loop over jets
  const float jet_min_pt = 15;
  value_jet_n = 0;
  std::vector<CaloJet> selectedJets;
  for (auto it = jets->begin(); it != jets->end(); it++) {
    if (it->pt() > jet_min_pt) {
      value_jet_btag[value_jet_n] = btags->operator[](it - jets->begin()).second;
      value_jet_n++;
    }
  }
~~~
{: .source}

{: .output}

We can look at the b tag information in a dataset using ROOT. 

~~~
root root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/GluGluToHToTauTau.root 
TBrowser b;
~~~
{: .source}

{: .output}

Double click Jet_btag to see the distribution of discriminator values.


The b jet is tagged if the discriminator value exceeds some threshold which is called the cut value.
The BTV POG defines operating working points for the b-tagging algorithms to define cuts in the discriminator and 
are based on a misidentification probability for light-flavor jets in ttbar events of:

* ~10% (Loose)
* ~1% (Medium)
* ~0.1% (Tight)

The discriminator cuts from working points defined by the BTV POG for this dataset are (https://twiki.cern.ch/twiki/bin/view/CMSPublic/BtagRecommendation2011OpenData):

CSVL (Loose) - 0.244
CSVM (Medium) - 0.679
CSVT (Tight) - 0.898

In the following code, input the discriminator cut for the medium WP:

~~~
for(unsigned int ijet=0; ijet < Jet_pt->size(); ijet++){

        if(Jet_pt->at(ijet) < jetPtCut || fabs(Jet_eta->at(ijet)) > jetEtaCut) continue;
        bool istagged =  Jet_btag > <PUT YOUR B-TAG DISCR CUT HERE>; 
        ...
}
~~~
{: .source}


Scale factors (ADD MORE):

 * [Explanation](https://twiki.cern.ch/twiki/bin/view/CMSPublic/BtagRecommendation2011OpenData#Data_MC_Scale_Factors)
 * [Data file for the CSV algoritm](https://twiki.cern.ch/twiki/pub/CMSPublic/BtagRecommendation2011OpenData/CSV.csv)
 * [Examples of application methods](https://twiki.cern.ch/twiki/bin/view/CMSPublic/BtagRecommendation2011OpenData#Methods_to_Apply_b_Tagging_Effic)
---
title: "CMS Jets and MET"
teaching: 15
questions:
- "How are jets and missing transverse energy accessed in CMS files?"
objectives:
- "Identify jet and MET code collections in AOD files"
- "Understand typical features of jet/MET objects"
- "Practice accessing jet quantities"
keypoints:
- "A" 
- "B"
- "C"
- "D"
---

Jets are spatially-grouped collections of long-lived particles that are produced when a quark or gluon hadronizes. The kinetmatic properties of
jets resemble that of the initial partons that produced them. In the CMS language, jets are made up of many particles, with the
following predictable energy composition:

*   ~65% charged hadrons
*   ~25% photons (from neutral pions)
*   ~10% neutral hadrons

Jets are very messy! Hadronization and the subsequent decays of unstable hadrons can produce 100s of particles near each other in the CMS detector.
Hence these particles are rarely analyzed individually. How can we determine which particle candidates should be included in each jet?

## Clustering

Jets can be clustered using a variety of different inputs from the CMS detector. "CaloJets" use only calorimeter energy deposits. "GenJets" use generated
particles from a simulation. But by far the most common are "PFJets", from **particle flow candidates**.

The result of the CMS Particle Flow algorithm is a list of particle candidates that account for all inner-tracker and muon tracks and all above-threshold
energy deposits in the calorimeters. These particles are formed into jets using a "clustering algorithm". The most common algorithm used by CMS is the
"anti-kt" algorithm, which is abbreviated "AK". It iterates over particle pairs and finds the two (*i* and *j*) that are the closest in some distance
measure and determines whether to combine them:

<a href="https://www.codecogs.com/eqnedit.php?latex=d_{ij}&space;=&space;min(p^{-2}_{T,i},p^{-2}_{T,j})\Delta&space;R^2_{ij}/R^2" target="_blank"><img src="https://latex.codecogs.com/svg.latex?d_{ij}&space;=&space;min(p^{-2}_{T,i},p^{-2}_{T,j})\Delta&space;R^2_{ij}/R^2" title="d_{ij} = min(p^{-2}_{T,i},p^{-2}_{T,j})\Delta R^2_{ij}/R^2" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\text{Combine&space;when&space;}&space;d_{ij}&space;<&space;p^{-2}_{T,i}\text{&space;;&space;stop&space;when&space;}&space;d_{ij}&space;>&space;p^{-2}_{T,i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\text{Combine&space;when&space;}&space;d_{ij}&space;<&space;p^{-2}_{T,i}\text{&space;;&space;stop&space;when&space;}&space;d_{ij}&space;>&space;p^{-2}_{T,i}" title="\text{Combine when } d_{ij} < p^{-2}_{T,i}\text{ ; stop when } d_{ij} > p^{-2}_{T,i}" /></a>

<img src="{{ page.root }}/clustering.png" alt="Clustering sequence" />

The momentum power (-2) used by the anti-kt algorithm means that higher-momentum particles are clustered first. This leads to jets with a round shape that
tend to be centered on the hardest particle. In CMS software this clustering is implemented using the [[fastjet][www.fastjet.fr]] package. 

<img src="{{ page.root }}/antikt.png" alt="anti-kt example" />

## Pileup

Inevitably, the list of particle flow candidates contains particles that did not originate from the primary interaction point. CMS experiences multiple
simultaneous collisions, called "pileup", during each "bunch crossing" of the LHC, so particles from multiple collisions coexist in the detector.
There are various methods to remove their contributions from jets:

 * Charged hadron subtraction (CHS): all charged hadron candidates are associated with a track. If the track is not associated with the primary vertex, that
 charged hadron can be removed from the list. CHS is limited to the region of the detector covered by the inner tracker. The pileup contribution to
 neutral hadrons has to be removed mathematically -- more in episode 3!
 * PileUp Per Particle Identification (PUPPI): CHS is applied, and then all remaining particles are weighted based on their likelihood of arising from
 pileup. This method is more stable and performant in high pileup scenarios such as the upcoming HL-LHC era.

## Accessing jets in CMS software

Jets software classes have the same basic 4-vector methods as the objects discussed in the previous lesson:

~~~
Handle<CaloJetCollection> jets;
iEvent.getByLabel(InputTag("ak5CaloJets"), jets);

for (auto it = jets->begin(); it != jets->end(); it++) {

    value_jet_pt[value_jet_n] = it->pt();
    value_jet_eta[value_jet_n] = it->eta();
    value_jet_phi[value_jet_n] = it->phi();
    value_jet_mass[value_jet_n] = it->mass();

}
~~~
{: .source}

Talk about the energy fractions, and how IDs are formed (FIND PUBLIC LINK!)

CHALLENGE: implement a "loose ID" boolean variable in teh code

## MET

Talk about MET

~~~
Handle<PFMETCollection> met;
iEvent.getByLabel(InputTag("pfMet"), met);

value_met_pt = met->begin()->pt();
value_met_phi = met->begin()->phi();
value_met_sumet = met->begin()->sumEt();

value_met_significance = met->begin()->significance();
auto cov = met->begin()->getSignificanceMatrix();
value_met_covxx = cov[0][0];
value_met_covxy = cov[0][1];
value_met_covyy = cov[1][1];

~~~
{: .source}

{% include links.md %}


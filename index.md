---
layout: lesson
root: .  # Is the only page that doesn't follow the pattern /:path/index.html
permalink: index.html  # Is the only page that doesn't follow the pattern /:path/index.html
---
Jets and missing transverse energy (MET) are critical for CMS physics analyses. They are more 
complex than most of the objects we discussed in the previous lesson, because they are 
reconstructed using **multiple particle-flow candidates**. After all candidates have been built
from the tracks and energy deposits in CMS, they can be "clustered" using a variety of 
algorithms into composite objects called "jets". Missing transverse energy clusters, in a sense, 
all candidates in the entire detector: it is the negative vector sum of the momentum of all
candidates. 

In this lesson we will explore the basic utilities for jets and MET, how to identify jets that 
arise from interesting original particles such as bottom quarks, and how to correct jets and 
MET for differences between data and simulation. 

<!-- this is an html comment -->

{% comment %} This is a comment in Liquid {% endcomment %}

> ## Prerequisites
>
> Complete the pre-exercises and watch videos X and Y and Z. 
{: .prereq}

> ## Big questions
> 
> 1. How are jets and missing transverse energy accessed in CMS files?
> 2. How can I identify jets from b quarks?
> 3. How are data/simulation differences dealt with for jets and MET?
> 4. <TBA>?
>
{: .checklist}

{% include links.md %}

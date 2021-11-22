# computerGraphicsCW
# Computer Graphics Coursework Assignment Brief

## Introduction

This document provides an overview of the coursework assignment for
**COMS30064 Computer Graphics**. This coursework is the only assignment for this unit
and as such is worth 15 credit points.

The assignment will be released on Monday 22<sup>nd</sup> November (Week 9) at 9am and your work
must be submitted by Friday 10<sup>th</sup> December (Week 11) at 1pm at the very latest.
The intention is that you submit by 12 noon and keep the last hour as emergency reserve for
(e.g. technical problems). In case of problems with your submission, you must e-mail
coms-info@bristol.ac.uk _before_ the 1pm final deadline to avoid your work being counted as late.

You must submit the coursework on the Blackboard page for the assessment unit **COMS30064**
(NOT the teaching unit **COMS30020**). On the assessment unit page, go to “Assessment, Submission
and Feedback” and follow the instructions found there.

You are expected to work on both of your optional unit courseworks in the 3-week period
from Week 9 to Week 11 as if it were a working week in a regular job - that is 5 days a week
for no more than 8 hours a day. It is up to you how you distribute your time and workload
between the two units within those constraints.

You are strongly advised NOT to try and work excessive hours during the coursework period:
this is more likely to make your health worse than to make your marks better. If you need
further pastoral/mental health support, please talk to your personal tutor, a senior tutor,
or the university wellbeing service.

For this unit **COMS30064 Computer Graphics**, online office hours will be provided on Teams
on Monday afternoons (2pm-4pm). In week 9 this session will be a whole-cohort Q&A session
(to answer any questions regarding the brief that was released that morning). In weeks 10 and 11
this session will be an online (Teams) drop-in office hour for individuals to ask questions.
Note that this online session is for high-level Q&A and as such is NOT intended as a
"debugging service" - if you are having problems getting your code to run you should use the
usual Teams discussion forum.

This coursework is an individual piece of work. You should therefore not collaborate with
other students - if you do, this will be considered as an academic offence and will be dealt
with in accordance with university plagiarism processes.

## Assignment Overview

The aim of this unit has been to provide you with an understanding
of a variety of different approaches to 3D rendering. Pre-recorded lecture segments
and narrated animations provided the main theoretical content of this unit
and the tasks from the weekly workbooks have allowed you to practically explore
the key concepts presented.

The assignment for this unit provides an opportunity to demonstrate your
understanding of these concepts, as well as a chance to exercise the codebase
that you have been accumulating over the course of the unit.

Many of the approaches to rendering that we have explored can only be fully
appreciated by movement of the camera through the model (or elements of the model itself).
For this reason, your task is to create a short (approximately 20 second)
animated "ident" video. For more insight into the concept of an "ident", you should view this
[compilation of BBC2 idents](https://www.youtube.com/watch?v=nPnXQMlq2iM).
Note that there are a wide variety of examples presented - some are
more relevant to this unit than others !

To provide a focus and purpose for your animation, you should create an ident
for promotional use by the
[MVB hackspace](https://www.facebook.com/groups/363509640788223/)

We have provided a texture-mapped hackspace logo which you may
like to use to help you on your way. This however is not mandatory and you
are free to use any alternative models you see fit.

## Assignment Details

You should refine and extend the codebase that you have developed during the
completion of the weekly workbooks in order to build a capable C++ renderer.
If you wish to attain a high mark, you should attempt to implement some of
the approaches introduced in the "advanced topics" playlist.

It is up to you to select which approaches you implement, however you should
remember that there is a performance trade-off at work. The more complex
the features you implement, the longer will be the rendering time for your
animation. You will need to find an appropriate compromise between:
renderer sophistication, algorithmic performance, animation complexity
animation aesthetics and rendering time.

You should use your code to render out individual frames of an animation (as
images) and then sequence these together into an mp4 using an external video
manipulation tool (guidance will be provided on suitable tools for this purpose).

You MUST ensure that your animation clearly demonstrates all of the rendering
approaches that you have implemented. If we can't see the feature operating correctly,
then we can't award marks for it. Note that you may find some features of your
renderer can only be fully demonstrated using certain modes of rendering
(wireframe/rasterising/raytracing) so you may need to switching between these during
your animation.

## Marking criteria

Your work will be assess on the extent to which you achieve the following objectives:

- Successful implementation of different approaches to 3D rendering
(wireframe, rasterising, raytracing)

- Successful implementation of a range of different approaches to lighting and shading
(ambient/diffuse/specular, hard/soft shadows, Phong/Gouraud shading)

- Effective animation of camera position, camera orientation and model elements

- Selection and implementation of advanced rendering features, such as:
    - "Exotic" materials (e.g. glass, mirrors, metal)
    - Surface mapping (texture maps, bump maps, environment maps)
    - Advanced lighting approaches (e.g. photon maps)

An indicative marking guide has been provided to illustrate the
the grades you might expect to receive for a given set of features and capabilities.

## Submission

You should submit a single 480p (640×480) MP4 video file containing your ident.
This MUST be the video file itself (i.e. not the URL to an online video).
You should make an effort to compress your video in order to reduce the file size.

You should also submit your code so that we may verify the correct implementation of
features. Include a Make file to build your rendered and ensure that your code compiles
on the MVB lab machines (this is where your code will be tested and marked !).
As with any programming assignment, you code will be checked against other
submissions using an automated plagiarism checker.

Upon submission, you must fill out
[this questionnaire](https://forms.office.com/Pages/ResponsePage.aspx?id=MH_ksn3NTkql2rGM8aQVG37h-tXTP5NGqEknVplKlalUNjJKREFLNTJIWEpMUFc5SkY5OEpQSDZKQy4u)
indicating which approaches you have used to implement specific rendering features. The aim of this
is to help optimise the marking process (so that we don't spend time searching for a
feature or approach that you haven't actually implemented !)

Students should submit all required materials to the "Assessment, submission and
feedback" section of Blackboard. Note that it is essential this is done on the
Blackboard page related to the "With Coursework" variant of the unit !

## Academic Offences

Academic offences (including submission of work that is not your own, falsification
of data/evidence or the use of materials without appropriate referencing) are all
taken very seriously by the University. Suspected offences will be dealt with in
accordance with the University’s policies and procedures. If an academic offence
is suspected in your work, you will be asked to attend a plagiarism panel, where
you will be given the opportunity to defend your work. The panel are able to apply
a range of penalties, depending the severity of the offence. These include:
requirement to resubmit assignment; capping of grades to the pass mark (40%);
awarding a mark of zero for an element of assessment; awarding a mark of zero for the
entire unit.

# sPHENIX-TPC-Lineshape

This repo show a few example pulse shape from the 2019 sPHENIX TPC test beam. 

[Lineshape.C](Lineshape.C) give the fit function `SignalShape_PowerLawDoubleExp()` 
that is based on Semi-Gaussian pulse shaping function, `signal_core = pow(t/tau, n) * exp( - t/tau)`. 
The fitting routine is isolated from the test beam analysis code: https://github.com/sPHENIX-Collaboration/prototype/tree/master/offline/packages/tpc2019

A few example ADC wavelet is fit with parameters shown in the same folder. 

# The data

Here we just load one output file from run 300, when we put 120 GeV proton beam through the center of the TPC prototye that measures the beam position in 16 planes. The signal is digitized on SAMPAv4 ASIC digitized at 20MHz. Event 10 where a track through the detector is selected and its clusterss are fit to the above function.

![Test beam](https://github.com/sPHENIX-Collaboration/tutorials/raw/0a12b21ed5b7950106d29e79597eba869d9404d0/JupyterLab/images/2019-sPHENIX-testbeam.jpg)
![Test beam](https://github.com/sPHENIX-Collaboration/tutorials/raw/0a12b21ed5b7950106d29e79597eba869d9404d0/JupyterLab/images/Beam-and-TPC.png)


Run 300 is from the Third Position Scan, (Lower Transfer Gap Voltage, better drift field) with following settings:


### Setup

Cathode: 20000 V

Gap:  150 V

GEM:  380 V

### Run list

<table>
<tbody>
<tr class="odd">
<td><p>Postion (in)</p></td>
<td><p>Runs</p></td>
</tr>
<tr class="even">
<td><p>6</p></td>
<td><p>297, 298<br />
</p></td>
</tr>
<tr class="odd">
<td><p>10</p></td>
<td><p>299<br />
</p></td>
</tr>
<tr class="even">
<td><p>14</p></td>
<td><p>300<br />
</p></td>
</tr>
<tr class="odd">
<td><p>18<br />
</p></td>
<td><p>301<br />
</p></td>
</tr>
</tbody>
</table>

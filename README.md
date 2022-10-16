# Gillespie

The following reaction scheme models a situation in which a modified gene produces two mRNAs which stochastically transcribe different reporter proteins, such as GFP and CFP, independently of one another forming the basis of a dual reporter method. The system features feedback from the mRNA given by its production rate being a function $f(x_{1})$.

$\text{Transcription: } \text{ } x_{1} \xrightarrow[]{f(x_{1})} x_{1} + 1$

$\text{Translation of } x_{2} \text{: } \text{ } x_{2} \xrightarrow[]{\lambda_{2} x_{1}} x_{2} + 1$

$\text{Translation of } x_{3} \text{: } x_{3} \xrightarrow[]{\lambda_{2} x_{1}} x_{3} + 1$

$\text{mRNA degradation: } \text{ } x_{1} \xrightarrow[]{\beta_{1} x_{1}} x_{1} - 1$

$\text{Degradation of } x_{2} \text{: } x_{2} \xrightarrow[]{\beta_{2} x_{2}} x_{2} - 1$

$\text{Degradation of } x_{3} \text{: } x_{3} \xrightarrow[]{\beta_{2} x_{3}} x_{3} - 1$

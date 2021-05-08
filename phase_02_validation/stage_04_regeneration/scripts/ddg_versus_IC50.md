${\text{The Gibbs free energy change of inhibition reaction before and after point mutation on peptide is expressed as}}$
$$
\Delta G(bind)=G_{comp}-G_{pro}-G_{pep}\\
\Delta G^{'}(bind)=G_{comp}^{'}-G_{pro}-G_{pep}^{'}
$$
${\text{We define}}~ddG~\text{to describe the effect of point mutation on binding affinity}$
$$
ddG=\Delta G^{'}(bind)-\Delta G(bind)
$$
$\text{The relationship between Gibbs free energy change and equilibrium constant is expressed as}$
$$
\Delta G(bind)=-RTln(K_c)
$$
${\text{When the inhibition reaches equilibrium, we have}}$
$$
K_c=\frac{c_{comp}}{c_{pro}c_{pep}}
$$
${\text{According to the definition of half maximal inhibitory concentration, we assume that}}$
$$
c_{comp}=c_{pro}~\text{when}~c_{pep}=IC_{50}
$$
${\text{Thus we have}}$
$$
ddG=RTln(\frac{c_{pep}^{'}}{c_{pep}})
$$
$\text{Substitute}~G~\text{with the value of score in Rosetta,}~t\text{, and we expect}$
$$
[(t_{comp}^{'}-t_{pep}^{'})-(t_{comp}-t_{pep})]\propto ln(\frac{c_{pep}^{'}}{c_{pep}})
$$
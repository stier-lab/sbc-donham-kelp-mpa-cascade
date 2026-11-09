"""Independent R/Python conversion agreement at representative sizes."""
from pathlib import Path
import subprocess,os,sys
import pandas as pd,numpy as np
A=Path(__file__).resolve().parents[2];D=A.parent/'sbc-donham-kelp-mpa-data';O=A/'outputs/integrated_review_20261103';sys.path.insert(0,str(A/'code/integrated_review'));from common import mass
code='''Sys.setenv(DONHAM_KELP_CALIBRATION_CSV="%s")
source("code/R/00c_analysis_constants.R")
source("%s")
x<-c(1,10,25,50,100,200)
write.csv(data.frame(size=x,lobster=bio_lobster(x),red=bio_redurch(x),purple=bio_purpurch(x),sheephead=bio_sheephead(x)),"%s",row.names=FALSE)
'''%(A/'data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv',D/'code/R/01_utils.R',O/'r_mass_verification.csv')
subprocess.run(['Rscript','-e',code],check=True,cwd=D)
r=pd.read_csv(O/'r_mass_verification.csv');rows=[]
for col,tx in [('lobster','P. interruptus'),('red','M. franciscanus'),('purple','S. purpuratus'),('sheephead','S. pulcher')]:
 expected=mass(tx,r['size']);assert np.allclose(expected,r[col],rtol=1e-12);rows.append(dict(Taxa=tx,largest_absolute_difference=float(np.max(np.abs(expected-r[col]))),passed=True))
pd.DataFrame(rows).to_csv(O/'conversion_verification.csv',index=False);print('R/Python conversions agree for all four animal taxa')

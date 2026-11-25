# Reference verification report

Course: *Longitudinal Data Analysis: Concepts, Methods and Practical Applications Using R*
Checked: 26 October 2026. Sources: Crossref metadata (api.crossref.org, which is what doi.org resolves to), PubMed E-utilities, CRAN package pages and CITATION files, and publisher, author and regulator pages.
BibTeX file: `references.bib` (only verified entries; citation key shown in brackets).

## Special question: "Longitudinal Data Analysis, Molenberghs and Verbeke 2017"

**Most likely referent:** Molenberghs, G. & Verbeke, G. (2017). *Longitudinal Data Analysis.* In *Wiley StatsRef: Statistics Reference Online*, pp. 1–28. John Wiley & Sons. Published online 15 May 2017. doi:10.1002/9781118445112.stat00291.pub2. Key: `[molenberghs2017]`.

- This is the only item with this exact title, these two authors in this order and the year 2017. It is an updated encyclopedia article. The earlier version of the same StatsRef entry (`stat00291`, 2014) was written by Burton Singer, and the 2017 update (`.pub2`) is by Molenberghs and Verbeke.
  Evidence: https://api.crossref.org/works/10.1002/9781118445112.stat00291.pub2 ; https://doi.org/10.1002/9781118445112.stat00291.pub2 ; earlier version: https://api.crossref.org/works/10.1002/9781118445112.stat00291
- The other candidates do not match 2017:
  - Fitzmaurice, Davidian, Verbeke & Molenberghs (eds), *Longitudinal Data Analysis*, Chapman & Hall/CRC Handbooks of Modern Statistical Methods. Crossref and the publisher give **2008** (doi:10.1201/9781420011579, ISBN 978-1-58488-658-7). Some library catalogues say 2009. I found **no 2017 or later edition or reprint**, and no paperback ISBN. Evidence: https://api.crossref.org/works/10.1201/9781420011579 ; https://documentserver.uhasselt.be/handle/1942/9642
  - Verbeke & Molenberghs, *Linear Mixed Models for Longitudinal Data*, Springer **2000** (doi:10.1007/b98969). A **2009 softcover reprint** does exist (ISBN 978-1-4419-0299-3; eBook 978-1-4419-0300-6; doi:10.1007/978-1-4419-0300-6). Crossref lists Molenberghs first for the reprint, but that looks like a metadata error: catalogues and UHasselt give Verbeke & Molenberghs. Evidence: https://api.crossref.org/works/10.1007/978-1-4419-0300-6 ; https://katalog.muni.cz/Record/MUB01006381482/Details ; https://documentserver.uhasselt.be/handle/1942/10776
  - Molenberghs & Verbeke, *Models for Discrete Longitudinal Data*, Springer **2005**. No 2017 edition was found.
  - Other 2017 papers by the pair: Verbeke & Molenberghs (2017) "Modeling Through Latent Variables", *Annual Review of Statistics and Its Application* 4:267–282 (doi:10.1146/annurev-statistics-060116-054017). Verbeke, Fieuws, Molenberghs & Davidian (2017), *SMMR* 26(1):112, a response to a letter about their 2014 review. Neither is titled "Longitudinal Data Analysis".
- **Recommendation:** Cite `molenberghs2017` if the course owner meant a 2017 source. If they meant "the Molenberghs/Verbeke longitudinal book", they most likely mean Verbeke & Molenberghs (2000; 2009 reprint) or the 2008 Fitzmaurice et al. handbook, and you should ask them which. All three are in the .bib file.

## Item-by-item results

Status: **V** = verified as given; **C** = corrected (change listed); **N** = note or caveat.

| # | Item | Key | Status | Corrections / notes | Evidence |
|---|------|-----|--------|---------------------|----------|
| 1 | Laird & Ware 1982 | laird1982 | V | doi:10.2307/2529876 | https://api.crossref.org/works/10.2307/2529876 ; PubMed 7168798 |
| 2 | Liang & Zeger 1986 | liang1986 | V | doi:10.1093/biomet/73.1.13 | https://api.crossref.org/works/10.1093/biomet/73.1.13 |
| 3 | Zeger, Liang & Albert 1988 | zeger1988 | V | doi:10.2307/2531734 | https://api.crossref.org/works/10.2307/2531734 ; PubMed 3233245 |
| 4 | Diggle et al. 2002, 2nd ed. | diggle2002 | V | Oxford Statistical Science Series 25; ISBN 9780198524847; doi:10.1093/oso/9780198524847.001.0001 | https://api.crossref.org/works/10.1093/oso/9780198524847.001.0001 ; https://eprints.lancs.ac.uk/id/eprint/51931/ |
| 5 | Fitzmaurice, Laird & Ware 2011, 2nd ed. | fitzmaurice2011 | V | ISBN 9780470380277; doi:10.1002/9781119513469 | https://api.crossref.org/works/10.1002/9781119513469 |
| 6 | Singer & Willett 2003 | singer2003 | V (N) | Full title adds the subtitle *Modeling Change and Event Occurrence*; doi:10.1093/acprof:oso/9780195152968.001.0001 | https://api.crossref.org/works/10.1093/acprof:oso/9780195152968.001.0001 |
| 7 | Hedeker & Gibbons 2006 | hedeker2006 | V | doi:10.1002/0470036486 | https://api.crossref.org/works/10.1002/0470036486 |
| 8 | Pinheiro & Bates 2000 | pinheiro2000 | V | Series Statistics and Computing; doi:10.1007/b98882 (2009 reprint: doi:10.1007/978-1-4419-0318-1) | https://api.crossref.org/works/10.1007/b98882 ; https://api.crossref.org/works/10.1007/978-1-4419-0318-1 |
| 9 | Bates et al. 2015 lme4 | bates2015 | V | pp. 1–48; doi:10.18637/jss.v067.i01 | https://api.crossref.org/works/10.18637/jss.v067.i01 ; https://cran.r-project.org/web/packages/lme4/citation.html |
| 10 | Kuznetsova et al. 2017 lmerTest | kuznetsova2017 | V | Title "lmerTest Package: Tests in Linear Mixed Effects Models", pp. 1–26; doi:10.18637/jss.v082.i13 | https://api.crossref.org/works/10.18637/jss.v082.i13 ; https://www.jstatsoft.org/v82/i13 |
| 11 | Brooks et al. 2017 glmmTMB | brooks2017 | V | Full title "glmmTMB Balances Speed and Flexibility Among Packages for Zero-inflated Generalized Linear Mixed Modeling"; 9 authors; doi:10.32614/RJ-2017-066 | https://api.crossref.org/works/10.32614/RJ-2017-066 |
| 12 | geepack (JSS 15(2)) | **halekoh2006** | **C** | **First author is Halekoh U, not Højsgaard S.** Correct order: Halekoh U, Højsgaard S, Yan J (2006), pp. 1–11, doi:10.18637/jss.v015.i02. (RePEc dates the issue Dec 2005; Crossref and CRAN give 2006.) | https://api.crossref.org/works/10.18637/jss.v015.i02 ; https://cran.r-project.org/web/packages/geepack/citation.html |
| 13 | van Buuren & Groothuis-Oudshoorn 2011 mice | vanbuuren2011 | V | pp. 1–67; doi:10.18637/jss.v045.i03 | https://api.crossref.org/works/10.18637/jss.v045.i03 ; https://www.jstatsoft.org/v45/i03 |
| 14 | van Buuren 2018, 2nd ed. | vanbuuren2018 | V | Chapman & Hall/CRC; doi:10.1201/9780429492259 | https://api.crossref.org/works/10.1201/9780429492259 |
| 15 | Little & Rubin 2019, 3rd ed. | little2019 | V | doi:10.1002/9781119482260 | https://api.crossref.org/works/10.1002/9781119482260 |
| 16 | Rubin 1976 | rubin1976 | V | doi:10.1093/biomet/63.3.581 | https://api.crossref.org/works/10.1093/biomet/63.3.581 |
| 17 | Kenward & Roger 1997 | kenward1997 | V | doi:10.2307/2533558 | https://api.crossref.org/works/10.2307/2533558 ; PubMed 9333350 |
| 18 | Mancl & DeRouen 2001 | mancl2001 | V | doi:10.1111/j.0006-341X.2001.00126.x | https://api.crossref.org/works/10.1111/j.0006-341X.2001.00126.x |
| 19 | Robins, Rotnitzky & Zhao 1995 | robins1995 | V | doi:10.1080/01621459.1995.10476493 | https://api.crossref.org/works/10.1080/01621459.1995.10476493 |
| 20 | Wood 2017, 2nd ed. | wood2017 | V | doi:10.1201/9781315370279 | https://api.crossref.org/works/10.1201/9781315370279 |
| 21 | Rizopoulos 2012 | rizopoulos2012 | V | doi:10.1201/b12208 | https://api.crossref.org/works/10.1201/b12208 |
| 22 | Curran & Bauer 2011 | curran2011 | V | doi:10.1146/annurev.psych.093008.100356 | https://api.crossref.org/works/10.1146/annurev.psych.093008.100356 |
| 23 | Enders & Tofighi 2007 | enders2007 | V (N) | Full title adds the subtitle "A new look at an old issue"; doi:10.1037/1082-989X.12.2.121 | https://api.crossref.org/works/10.1037/1082-989X.12.2.121 |
| 24 | Hernán & Robins 2020 | hernan2020 | V (N) | The authors' recommended citation is "Hernán MA, Robins JM (2020). Causal Inference: What If. Boca Raton: Chapman & Hall/CRC". The book is online and has no DOI. The authors' page says the print edition has not been released yet. The CRC page (ISBN 9781420076165) shows a 1st edition of 312 pages with no confirmed print date. | https://miguelhernan.org/whatifbook ; https://www.routledge.com/Causal-Inference-What-If/Hernan-Robins/p/book/9781420076165 |
| 25 | Molenberghs & Kenward 2007 | molenberghs2007 | V | doi:10.1002/9780470510445 | https://api.crossref.org/works/10.1002/9780470510445 |
| 26 | National Research Council 2010 | nrc2010 | V | doi:10.17226/12955; ISBN 9780309158145 | https://api.crossref.org/works/10.17226/12955 |
| 27 | Self & Liang 1987 | self1987 | V | doi:10.1080/01621459.1987.10478472 | https://api.crossref.org/works/10.1080/01621459.1987.10478472 |
| 28 | Stram & Lee 1994 | stram1994 | V | doi:10.2307/2533455 | https://api.crossref.org/works/10.2307/2533455 ; PubMed 7786999 |
| 29 | Barr et al. 2013 | barr2013 | V | doi:10.1016/j.jml.2012.11.001 | https://api.crossref.org/works/10.1016/j.jml.2012.11.001 |
| 30 | Matuschek et al. 2017 | matuschek2017 | V | doi:10.1016/j.jml.2017.01.001 | https://api.crossref.org/works/10.1016/j.jml.2017.01.001 |
| 31 | Senn 2006 | senn2006 | V | doi:10.1002/sim.2682 | https://api.crossref.org/works/10.1002/sim.2682 |
| 32 | Vickers & Altman 2001 | vickers2001 | V | Issue 7321; doi:10.1136/bmj.323.7321.1123 | https://api.crossref.org/works/10.1136/bmj.323.7321.1123 |
| 33 | Carpenter & Kenward 2013 | carpenter2013 | V | doi:10.1002/9781119942283 | https://api.crossref.org/works/10.1002/9781119942283 |
| 34 | Pullenayegum & Lim 2016 | pullenayegum2016 | V | doi:10.1177/0962280214536537 | https://api.crossref.org/works/10.1177/0962280214536537 |
| 35 | Morris, White & Crowther 2019 | morris2019 | V | doi:10.1002/sim.8086 | https://api.crossref.org/works/10.1002/sim.8086 |
| 36 | Arnold et al. 2011 | arnold2011 | V | Article 94; doi:10.1186/1471-2288-11-94 | https://api.crossref.org/works/10.1186/1471-2288-11-94 |
| 37 | Gamble et al. 2017 | gamble2017 | V | 15 authors; doi:10.1001/jama.2017.18556 | https://api.crossref.org/works/10.1001/jama.2017.18556 ; PubMed 29260229 |
| 38 | ICH E9(R1) 2019 | ich2019 | V | Step 4 adopted 20 Nov 2019 (document dated 17 Nov 2019); EMA/CHMP/ICH/436221/2017; EMA effective date 30 Jul 2020 | https://www.ema.europa.eu/documents/scientific-guideline/ich-e9-r1-addendum-estimands-sensitivity-analysis-clinical-trials-guideline-statistical-principles_en.pdf ; https://www.pmda.go.jp/files/000269157.pdf |
| 39 | Lindstrom & Bates 1990 | lindstrom1990 | V | doi:10.2307/2532087 | https://api.crossref.org/works/10.2307/2532087 ; PubMed 2242409 |
| 40 | Daniels & Hogan 2008 | daniels2008 | V | doi:10.1201/9781420011180 | https://api.crossref.org/works/10.1201/9781420011180 |
| 41 | Gelman & Hill 2007 | gelman2007 | V (N) | Crossref gives a print date of 18 Dec 2006, but the copyright page says first published 2007, so 2007 is kept. Series: Analytical Methods for Social Research; doi:10.1017/CBO9780511790942 | https://api.crossref.org/works/10.1017/CBO9780511790942 ; https://assets.cambridge.org/97805218/67061/copyright/9780521867061_copyright_info.pdf |
| 42 | Wickham 2014 Tidy data | wickham2014 | V (N) | Vol. 59, issue 10 and the DOI 10.18637/jss.v059.i10 are verified. **The page range could not be verified** (JSS site unreachable), so the .bib entry has no pages field. | https://api.crossref.org/works/10.18637/jss.v059.i10 ; https://ideas.repec.org/a/jss/jstsof/v059i10.html |
| 43 | Greenhouse & Geisser 1959 | greenhouse1959 | V | doi:10.1007/BF02289823 | https://api.crossref.org/works/10.1007/BF02289823 |
| 44 | Huynh & Feldt 1976 | huynh1976 | V | doi:10.3102/10769986001001069 | https://api.crossref.org/works/10.3102/10769986001001069 |
| 45 | Heinze & Schemper 2002 | heinze2002 | V | doi:10.1002/sim.1047 | https://api.crossref.org/works/10.1002/sim.1047 |
| 46 | DHARMa | hartig2026 | V | Current CRAN version **0.5.0** (published 2026-06-01). There is no CITATION file, so R's default citation applies: Hartig F (2026). DHARMa: Residual Diagnostics for Hierarchical (Multi-Level / Mixed) Regression Models. R package version 0.5.0. | https://cran.r-project.org/web/packages/DHARMa/index.html |
| 47 | emmeans | lenth2026 | **C** | Version **2.0.4** (2026-07-15). The authors listed as "aut" are now **Lenth RV and Piaskowski J**, so the default citation is "Lenth R, Piaskowski J (2026)", not Lenth alone. There is no CITATION file. | https://cran.r-project.org/web/packages/emmeans/index.html |
| 48 | ordinal | christensen2026 | V | Package CITATION: Christensen RHB (2026). ordinal---Regression Models for Ordinal Data. R package version 2026.7-26 | https://cran.r-project.org/web/packages/ordinal/citation.html |
| 49 | JMbayes2 | rizopoulos2026 | **C** | **Author order on CRAN is Rizopoulos D, Miranda-Afonso P, Papageorgiou G** (the request had Papageorgiou second). Version 0.6-5 (2026-10-15). There is no CITATION file. | https://cran.r-project.org/web/packages/JMbayes2/index.html |
| 50 | renv | ushey2026 | V | Ushey K, Wickham H; version 1.3.1 (2026-10-06); copyright holder and funder Posit Software, PBC | https://cran.r-project.org/web/packages/renv/index.html |
| 51 | Neuhaus, Kalbfleisch & Hauck 1991 | neuhaus1991 | V (N) | Full title ends "...for Analyzing Correlated Binary Data"; doi:10.2307/1403572. Crossref gives only the first page (25). The last page (35) comes from secondary citations. | https://api.crossref.org/works/10.2307/1403572 ; https://www.stata.com/meeting/3uk/slides/ |
| 52 | Verbeke & Molenberghs 2000 | verbeke2000 | V | doi:10.1007/b98969; 2009 softcover reprint ISBN 978-1-4419-0299-3 (doi:10.1007/978-1-4419-0300-6) | https://api.crossref.org/works/10.1007/b98969 ; https://documentserver.uhasselt.be/handle/1942/10776 |
| 53 | Molenberghs & Verbeke 2005 | molenberghs2005 | V | Springer Series in Statistics; doi:10.1007/0-387-28980-1 | https://api.crossref.org/works/10.1007/0-387-28980-1 |
| 54 | Fitzmaurice et al. (eds) handbook | fitzmaurice2008 | **C** (year added) | **Year 2008** (Crossref online date 11 Aug 2008; some catalogues say 2009). Series: Chapman & Hall/CRC Handbooks of Modern Statistical Methods; ISBN 9781584886587; doi:10.1201/9781420011579 | https://api.crossref.org/works/10.1201/9781420011579 |
| 55 | Twisk et al. 2018 | twisk2018 | V | Authors: Twisk J, Bosman L, Hoekstra T, Rijnhart J, Welten M, Heymans M; doi:10.1016/j.conctc.2018.03.008 | https://api.crossref.org/works/10.1016/j.conctc.2018.03.008 |
| 56 | Mallinckrodt et al. 2008 | mallinckrodt2008 | **C** (issue added) | **Issue 4** (vol 42(4):303–319). Authors: Mallinckrodt CH, Lane PW, Schnell D, Peng Y, Mancuso JP; doi:10.1177/009286150804200402 | https://api.crossref.org/works/10.1177/009286150804200402 |

## Items that could not be fully verified
- Wickham (2014) page range: omitted from the .bib entry. Everything else in that entry is verified.
- Neuhaus et al. (1991) last page (35): confirmed only through secondary sources. The volume, issue, first page and DOI are confirmed by Crossref.
- Hernán & Robins: no print publication date could be confirmed. The entry uses the authors' own citation (2020).
No item had to be left out of the .bib file.

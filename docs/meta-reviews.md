<!-- ------- Jan 2026 Cycle --------- -->
Paper Decision
Decisionby Program Chairs05 Apr 2026, 20:45 (modified: 06 Apr 2026, 19:31)Program Chairs, Area Chairs, Reviewers, AuthorsRevisions
Decision: Reject
Comment:
Your submission was recognized as a borderline case, but given the high volume of similarly recommended papers, it could not be accommodated in the final program.

Meta Review of Submission5910 by Area Chair jvJM
Meta Reviewby Area Chair jvJM26 Mar 2026, 12:19 (modified: 06 Apr 2026, 22:00)Area Chairs, Authors, Reviewers, Program ChairsRevisions
Metareview: 
This paper presents EngTrace, a benchmark for evaluating LLMs on physically grounded engineering reasoning. The benchmark is valuable. There are some necessary revisions that could be made for publication.

Confidence: 5: The area chair is absolutely certain
Rating: 5: Marginally below acceptance threshold
Recommendation: Borderline Findings
Presentation Mode: Poster


<!-- ------- May 2026 Cycle --------- -->
Meta Review of Submission2937 by Area Chair vT1J
Meta Reviewby Area Chair vT1J27 Jul 2026, 18:02 (modified: 03 Aug 2026, 19:17)Senior Area Chairs, Area Chairs, Authors, Reviewers Submitted, Program Chairs, Area Chair vT1J, Commitment ReadersRevisions
Metareview:
This paper introduces EngTrace, a symbolic benchmark for evaluating the engineering reasoning of large language models. It builds parameterized, contamination-resistant problem templates across several engineering branches and pairs each instance with a gold-standard reasoning trace, then proposes a two-stage evaluation framework that verifies intermediate reasoning steps rather than only final answers and applies it to a broad set of models. The work targets resources and evaluation for physically grounded reasoning.

Reviewers agree the symbolic-template design is a sound mechanism for contamination resistance and gold-standard traces, and the process-level evaluation with human-alignment validation is a reasonable and useful direction. The main concerns were a circularity risk because the judging panel overlaps with the evaluated models, insufficient statistical rigor with single-run scores and no variance or significance testing, under-supported causal claims about mathematical pretraining that confound specialization with scale, limited synthetic scope with no tool-augmented baselines, and presentation issues including a broken cross-reference. The rebuttal responded substantively, running the variance and bootstrap significance analysis, pointing to an existing threshold-sensitivity table, providing a same-family same-scale controlled comparison, and committing to a direct judge-exclusion ablation and broader coverage in the camera-ready version, though much of this remains deferred and no reviewer changed their score.

Overall, this is a useful benchmark whose core contribution is valuable, but several central claims still depend on additions promised for the camera-ready version and the reviews remain at the borderline with no score changes.

Summary Of Reasons To Publish:
The benchmark addresses a genuinely under-served need by evaluating physically grounded engineering reasoning at the level of intermediate steps rather than only final answers, using symbolic templates that resist contamination and yield gold-standard traces.
The evaluation covers a broad set of models and surfaces informative findings, including qualitatively distinct failure modes between frontier and open-weight models and a complexity cliff on harder problems.
The process-level metric is validated against a blind human-expert study rather than merely asserted, which is the appropriate kind of evidence for this class of contribution.
Summary Of Suggested Revisions:
Fold the rebuttal additions into the paper: report variance and significance testing for headline claims, and add the promised judge-exclusion ablation to directly address the panel-overlap circularity concern.
Strengthen the causal claim about mathematical pretraining with branch-level reporting for the math-specialized models, and promote key threshold-sensitivity and validation results from the appendices into the main text.
Clarify the framing so the synthetic scope and the physical-versus-linguistic diversity distinction are explicit, and complete a proofreading pass to fix the broken cross-reference and related typographic issues.
Overall Assessment: 2.5 = Borderline Findings
Reported Issues: No
Publication Ethics Policy Compliance: I did not use any generative AI tools for this review

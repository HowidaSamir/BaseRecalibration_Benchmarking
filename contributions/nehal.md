1. Teamleader: I have tried with my colleagues to assign tasks and make sure clear communication is established between team members. We had fruitful brainstorming to every step of the project starting from the dataset choice till the very last conclusion we could reach.
2. Going thoroughly through GATK Best Practice to make sure we choose the suitable best practice according to the type of the dataset.
https://software.broadinstitute.org/gatk/best-practices/workflow?id=11145
https://software.broadinstitute.org/gatk/best-practices/workflow?id=11146
https://software.broadinstitute.org/gatk/best-practices/workflow?id=11164
3. Investigating the feasibility of choosing E.coli to apply variant calling on and justified the reasons for my colleagues.
https://gatkforums.broadinstitute.org/gatk/discussion/10913/how-to-run-the-pathseq-pipeline
4. Searching for suitable dataset to start our project with. Learnt how to navigate through SRA and use its advance search options. The search options were explained to my colleagues to help them find datasets too.
5. Suggested a large group of different datasets to my colleagues.
6. Moderated my colleagues to share any educational resources they find online with each other.
7. Tried to help my colleague Eman working on the dataset she downloaded, however, hard disk space did not help.
8. Suggested the command  to count the FASTQ file reads(after searching): echo $(cat myfile.fastq|wc -l)/4|bc
9. Worked with my colleague Nada to index the VCF file of the whole human genome. Supplied the commands to change the name of the chromosome in the VCF file to be compatible with the reference file using the code:
grep "^#" Homo_sapiens.vcf > Homo_adjusted.vcf
perl -pe 's/^([^#])/chr\1/' Homo_sapiens.vcf >> Homo_adjusted.vcf 
10. Troubleshooted unparsed VCF files and shared the way with my colleagues on GitHub if we would like to remove a line:
    *To check for the line number:
    grep -nr "yourstring" yourfile.vcf
    *To delete the line:
    sed -i 'line#d' yourfile.vcf
    *To make sure line is deleted:
    grep -w "yourstring" yourfile.vcf
In the error message when it says "approximately", that doesn't have to be the exact line number and that's why we need to investigate for the line number first.

11. Writing the proposal introduction, aim, shared in materials, challenges, recommendations, and conclusion.
12. Fixed corrupted proposal file and supplying missing parts.
12. Extra tips for anyone working on a group project:
- Effective communication is the key for any project to stand and survive any challenges. By effective communication, I mean each group member has a solid idea of what other members are doing and ready to offer help as long as they can.
- Prior to starting to assign tasks to each member, especially when it comes to Bioinformatics, computational resources should be clearly clarified so that time and effort are saved.
- Computational resources include, but not limited to , Ubuntu operating system, Internet speed and accessbility, RAM, and memory disk space.
- Always keep draft for commands and every little detail about errors one faced and how they troubleshooted it..
- Metadata should be thoroughly examined and checked for suitability in pipeline.
- While troubleshooting, whether in a team or online, question should be stated clearly togeether with the faulty command and error received.
- Always keep a draft for everything. Things  get messed up and we are human and forget.


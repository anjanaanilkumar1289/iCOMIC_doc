## 6. Creating a custom pipeline

#### 6.1. Shell scripts to be written by the user
iCOMIC integrates most of the best practise tools for Whole Genome sequencing and RNA seq data analysis. The toolkit provides the user the complete freedom to choose any compatible combination of tools for analysis. iCOMIC permits a user to add a new tool as well. If you are a developer, adding a new tool is easy as iCOMIC relies on Snakemake where the codes are very readable. The primary thing to be done for integrating an additional tool is to create a rule file inside the directory `iCOMIC/rules`. Name the file as `[TOOL_NAME].smk`. The major thing to be taken care of is the name of input and output files. It should match the other tools of the corresponding analysis step.

![Figure: rule.smk](https://lh3.googleusercontent.com/Kt8r4GsAY0zpY6OMnyyTb8-EesQfQPhPYEdEjpQU05bSzT2cX4_YSmeI4BBP5tZZAk6D2S85iSPW1YYQooqYaAb1_jPALyYGYRK85oNoXNHjmN6HEDtS7OpveeMmCvKpYooCkhY5J6zh7VkZmRBw9MmEp0iQpOG0qwILZjFaQqJHcIgqkelm8VmiFu7y-NnQr02vTqdJQTQ1bRAD9_gaYx5YtqD9FTwtsUG12mn1g_WaVAns2ucVPPwUDf0dSXHcA7e4vJAgaJGmOe-UBW22BR_0LKU8dgKeZ-iZzqMym9l1P7JrSoWHipa3zidqcIXM-i4by7kmZX6sCqEbO-aiNXuRCX0U5YKd1T1JL14ZxxvkKlLPWMfqqSAFtFCAh8NdyfENAxSo9k492KTN-WveL7XlUzHdi5ah-tJTF9cSloQsXOLNta7S6pRbJyRQDuOA1Qk6YxfzCzObLLuyE_A8QAhJDzwkQ596SpXI4VJ1tIM91hhi1sGzBPajfBZrajTVWzEULreWPIrnkjwak7mYI8V95FOaba0BP_pey5GnDD1D4-UgGXpGOQfHb_u3Lk5U1X-ydRLg-9-Rklo9a0pfIeZvycUGviP0xNduke6TXOlAOUMPHx4N0aKFtOCXC6N1iI86ksgrTp20Dy51tc0kfZ_ZCCXTXcpAhQlldNNmVI-GbjqIQ9QTwGnrVXCcHz0=w893-h243-no?authuser=0)

This is how a rule file would look like. Parameters can be specified in the rule itself in the section `params`. If a snakemake wrapper is available for your choice of tool, that can be used, otherwise you need to write a shell command. List of snakemake wrappers are available in [Snakemake wrapper repository](https://snakemake-wrappers.readthedocs.io/en/stable/index.html). 

#### 6.2. How to run a custom pipeline

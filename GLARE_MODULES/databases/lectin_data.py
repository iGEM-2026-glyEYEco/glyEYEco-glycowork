
LECTIN_DATABASE = {

    "Gal-1": {
        "abbreviation":    ["Gal-1", "Galectin-1", "LGALS1"],
        "name":    "Galectin-1",
        "gene":    "LGALS1",
        "family":  "Galectin",
        "uniprot": "P09382",
        "primary": ["Gal(b1-4)GlcNAc", "Gal(b1-3)GalNAc"],
        "secondary":["Lac-NAc poly-N-acetyllactosamine", "Gal(b1-3)GlcNAc"],
        "negative": [],
        "tissues": ["T-cells", "smooth muscle", "endothelium", "thymus", "placenta"],
        "notes":   "Homodimeric proto-type galectin; pro-apoptotic on activated T-cells; "
                   "binds Core 1 and LacNAc termini. Key role in tumour immune evasion.",
        "doi": [
            "https://doi.org/10.1021/bi201121m",
            "https://doi.org/10.1093/glycob/cww002",
        ],
    },

    "Gal-2": {
        "abbreviation":    ["Gal-2", "Galectin-2", "LGALS2"],
        "name":    "Galectin-2",
        "gene":    "LGALS2",
        "family":  "Galectin",
        "uniprot": "P05162",
        "primary": ["Gal(b1-4)GlcNAc", "Lactose"],
        "secondary":["Gal(b1-3)GalNAc", "Gal(b1-3)GlcNAc"],
        "negative": [],
        "tissues": ["stomach — gastric mucosa", "small intestine"],
        "notes":   "Proto-type galectin; highest expression in stomach. "
                   "Involved in gastric mucosal defense; Helicobacter pylori modulates "
                   "gastric galectin expression. Less extensively studied than Gal-1/3.",
        "doi": [
            "https://doi.org/10.1093/glycob/cww002",
            "https://doi.org/10.1016/j.bbagen.2016.06.001",
        ],
    },

    "Gal-3": {
        "abbreviation":    ["Gal-3", "Galectin-3", "LGALS3", "Mac-2"],
        "name":    "Galectin-3",
        "gene":    "LGALS3",
        "family":  "Galectin",
        "uniprot": "P17931",
        "primary": ["Gal(b1-4)GlcNAc", "Gal(b1-3)GalNAc", "Neu5Ac(a2-3)Gal(b1-4)GlcNAc"],
        "secondary":["poly-N-acetyllactosamine", "Gal(b1-3)GlcNAc", "Gal(b1-4)Glc"],
        "negative": [],
        "tissues": ["macrophages", "monocytes", "cornea", "gut epithelium", "liver", "lung"],
        "notes":   "Chimera-type galectin; N-terminal aggregation domain + C-terminal CRD. "
                   "Binds MUC1 Core 1 in cancer. Elevated in fibrosis and ocular surface disease. "
                   "Serum marker in heart failure.",
        "doi": [
            "https://doi.org/10.1074/jbc.M213003200",
            "https://doi.org/10.1038/nrc1251",
            "https://doi.org/10.1016/j.celrep.2021.109819",
        ],
    },

    "Gal-4": {
        "abbreviation":    ["Gal-4", "Galectin-4", "LGALS4"],
        "name":    "Galectin-4",
        "gene":    "LGALS4",
        "family":  "Galectin",
        "uniprot": "P56470",
        "primary": ["Gal(b1-4)GlcNAc", "Gal(b1-3)GalNAc", "blood group A/B antigens"],
        "secondary":["sulfated glycosphingolipids", "3-O-sulfo-Gal"],
        "negative": [],
        "tissues": ["small intestine", "colon", "stomach epithelium"],
        "notes":   "Tandem-repeat galectin; two CRDs with distinct but overlapping specificities. "
                   "Stabilises lipid rafts in intestinal epithelium.",
        "doi": [
            "https://doi.org/10.1021/bi201121m",
            "https://doi.org/10.1093/nar/gkab918",
        ],
    },
    
    "Gal-7": {
        "abbreviation":    ["Gal-7", "Galectin-7", "LGALS7"],
        "name":    "Galectin-7",
        "gene":    "LGALS7",
        "family":  "Galectin",
        "uniprot": "P47929",
        "primary": ["Gal(b1-4)GlcNAc", "Gal(b1-3)GalNAc"],
        "secondary":["poly-N-acetyllactosamine"],
        "negative": [],
        "tissues": ["keratinocytes", "skin epidermis", "corneal epithelium"],
        "notes":   "p53-induced gene; regulates apoptosis in stratified squamous epithelium. "
                   "Expressed in cornea — relevant to ocular surface biology.",
        "doi": [
            "https://doi.org/10.1021/bi201121m",
        ],
    },
    "Gal-8": {
        "abbreviation":    ["Gal-8", "Galectin-8", "LGALS8"],
        "name":    "Galectin-8",
        "gene":    "LGALS8",
        "family":  "Galectin",
        "uniprot": "O00214",
        "primary": ["Neu5Ac(a2-3)Gal(b1-4)GlcNAc", "3-O-sulfo-Gal(b1-4)GlcNAc"],
        "secondary":["Gal(b1-4)GlcNAc", "poly-N-acetyllactosamine"],
        "negative": ["Neu5Ac(a2-6)Gal(b1-4)GlcNAc"],
        "tissues": ["liver", "kidney", "lung", "placenta", "T-cells"],
        "notes":   "Tandem-repeat galectin; N-CRD strongly prefers alpha2-3 sialylated and "
                   "3-O-sulfated LacNAc. C-CRD is broader. Activates integrins. "
                   "Role in autophagy of damaged lysosomes.",
        "doi": [
            "https://doi.org/10.1021/acschembio.1c00689",
            "https://doi.org/10.1074/jbc.M213003200",
            "https://www.sciencedirect.com/science/article/pii/S0021925820537587",
        ],
    },

    "Gal-9": {
        "abbreviation":    ["Gal-9", "Galectin-9", "LGALS9"],
        "name":    "Galectin-9",
        "gene":    "LGALS9",
        "family":  "Galectin",
        "uniprot": "O00182",
        "primary": ["Gal(b1-4)GlcNAc", "poly-N-acetyllactosamine"],
        "secondary":["Gal(b1-3)GalNAc", "Neu5Ac(a2-3)Gal"],
        "negative": [],
        "tissues": ["T-cells", "fibroblasts", "liver sinusoidal endothelium", "kidney"],
        "notes":   "Tandem-repeat galectin; ligand for Tim-3 (HAVCR2). Induces apoptosis "
                   "of Th1 and Th17 cells. Immune checkpoint target.",
        "doi": [
            "https://doi.org/10.1021/bi201121m",
            "https://doi.org/10.1093/glycob/cww002",
        ],
    },

    "Siglec-1": {
        "abbreviation":    ["Siglec-1", "CD169", "Sialoadhesin"],
        "name":    "Sialoadhesin / Siglec-1",
        "gene":    "SIGLEC1",
        "family":  "Siglec",
        "uniprot": "Q9BZZ2",
        "primary": ["Neu5Ac(a2-3)Gal(b1-4)GlcNAc", "Neu5Ac(a2-3)Gal(b1-3)GalNAc"],
        "secondary":["Neu5Gc(a2-3)Gal"],
        "negative": ["Neu5Ac(a2-6)"],
        "tissues": ["macrophages", "bone marrow sinusoidal macrophages", "spleen red pulp"],
        "notes":   "Exclusively alpha2-3 sialic acid. 17-domain ectodomain; mediates erythrocyte "
                   "and pathogen rosetting on macrophages.",
        "doi": [
            "https://doi.org/10.1073/pnas.2026102118",
            "https://doi.org/10.3390/cells10051260",
        ],
    },

    "Siglec-2": {
        "abbreviation":    ["Siglec-2", "CD22", "BL-CAM"],
        "name":    "Siglec-2 (CD22)",
        "gene":    "SIGLEC2",
        "family":  "Siglec",
        "uniprot": "P20273",
        "primary": ["Neu5Ac(a2-6)Gal(b1-4)GlcNAc", "Neu5Ac(a2-6)Gal(b1-4)Glc"],
        "secondary":["Neu5Ac(a2-6)GalNAc"],
        "negative": ["Neu5Ac(a2-3)Gal"],
        "tissues": ["B-cells (mature, follicular)", "B-cell lymphoma"],
        "notes":   "Inhibitory co-receptor on B-cells; dampens BCR signalling via SHP-1. "
                   "alpha2-6 sialylation requirement is strict — key difference from Siglec-9. "
                   "Target for anti-CD22 therapy (Inotuzumab, Epratuzumab) in B-ALL/lymphoma.",
        "doi": [
            "https://doi.org/10.1038/nri2202",
            "https://doi.org/10.1016/j.immuni.2006.05.011",
        ],
    },

    "Siglec-3": {
        "abbreviation":    ["Siglec-3", "CD33"],
        "name":    "CD33 / Siglec-3",
        "gene":    "CD33",
        "family":  "Siglec",
        "uniprot": "P20138",
        "primary": ["Neu5Ac(a2-3)Gal", "Neu5Ac(a2-6)Gal"],
        "secondary":["Neu5Gc"],
        "negative": [],
        "tissues": ["myeloid progenitors", "monocytes", "microglia (brain)"],
        "notes":   "Inhibitory ITIM receptor. Therapeutic target in AML (gemtuzumab). "
                   "Implicated in Alzheimer's disease via microglia.",
        "doi": [
            "https://doi.org/10.3390/cells10051260",
            "https://doi.org/10.1073/pnas.2026102118",
        ],
    },

    "Siglec-4": {
        "abbreviation":    ["Siglec-4", "MAG", "myelin-associated glycoprotein"],
        "name":    "Siglec-4 / MAG",
        "gene":    "MAG",
        "family":  "Siglec",
        "uniprot": "P20916",
        "primary": ["Neu5Ac(a2-3)Gal(b1-3)GalNAc", "3-sulfo-Gal",
                    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc"],
        "secondary":["GM1 ganglioside", "GD1a ganglioside", "GT1b ganglioside"],
        "negative": ["Neu5Ac(a2-6)"],
        "tissues": ["CNS myelin (Schwann cells / oligodendrocytes)", "peripheral nerve myelin"],
        "notes":   "Maintains axon-myelin stability; inhibits axon regeneration after injury. "
                   "Anti-MAG neuropathy is an autoimmune demyelinating disease driven by "
                   "antibodies to the 3-sulfo-glucuronyl (HNK-1) motif on MAG.",
        "doi": [
            "https://doi.org/10.1038/nrn2438",
            "https://doi.org/10.1016/j.neuron.2007.09.031",
        ],
    },

    "Siglec-7": {
        "abbreviation":    ["Siglec-7", "CD328", "p75/AIRM1"],
        "name":    "Siglec-7",
        "gene":    "SIGLEC7",
        "family":  "Siglec",
        "uniprot": "Q9Y286",
        "primary": ["Neu5Ac(a2-8)Neu5Ac", "GD3 ganglioside (Neu5Ac-a2-8-Neu5Ac-a2-3-Gal)"],
        "secondary":["Neu5Ac(a2-3)Gal(b1-4)GlcNAc"],
        "negative": ["Neu5Ac(a2-6)"],
        "tissues": ["NK cells", "monocytes", "CD8+ T-cells"],
        "notes":   "Inhibitory NK receptor; di-alpha2-8 sialic acid preference. Tumour cells "
                   "upregulate GD3 to evade NK killing.",
        "doi": [
            "https://doi.org/10.1073/pnas.2026102118",
            "https://doi.org/10.3390/cells10051260",
        ],
    },

    "Siglec-8": {
        "abbreviation":    ["Siglec-8"],
        "name":    "Siglec-8",
        "gene":    "SIGLEC8",
        "family":  "Siglec",
        "uniprot": "Q9NR61",
        "primary": ["Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
                    "6-sulfo-sialyl Lewis X"],
        "secondary":["Neu5Ac(a2-3)Gal(b1-4)GlcNAc"],
        "negative": [],
        "tissues": ["eosinophils", "mast cells", "basophils"],
        "notes":   "Highly restricted to eosinophils and mast cells. 6-sulfo sialyl Lewis X "
                   "is the preferred ligand. Therapeutic target in eosinophilic asthma.",
        "doi": [
            "https://doi.org/10.1073/pnas.2026102118",
            "https://doi.org/10.3390/cells10051260",
        ],
    },

    "Siglec-9": {
        "abbreviation":    ["Siglec-9", "CD329"],
        "name":    "Siglec-9",
        "gene":    "SIGLEC9",
        "family":  "Siglec",
        "uniprot": "Q9Y336",
        "primary": ["Neu5Ac(a2-3)Gal(b1-4)GlcNAc", "Neu5Ac(a2-6)Gal(b1-4)GlcNAc"],
        "secondary":["Neu5Ac(a2-3)Gal(b1-3)GalNAc"],
        "negative": [],
        "tissues": ["neutrophils", "monocytes", "NK cells", "dendritic cells"],
        "notes":   "Broadly expressed inhibitory Siglec; accepts both alpha2-3 and alpha2-6 "
                   "linkages. Immune checkpoint on neutrophils in tumour microenvironment.",
        "doi": [
            "https://doi.org/10.1073/pnas.2026102118",
            "https://doi.org/10.1016/j.sbi.2022.102412",
        ],
    },

    "Siglec-10": {
        "abbreviation":    ["Siglec-10", "CD170"],
        "name":    "Siglec-10",
        "gene":    "SIGLEC10",
        "family":  "Siglec",
        "uniprot": "Q9Y336",
        "primary": ["Neu5Ac(a2-3)Gal(b1-4)GlcNAc", "Neu5Ac(a2-6)Gal(b1-4)GlcNAc"],
        "secondary":["Neu5Gc"],
        "negative": [],
        "tissues": ["NK cells", "B-cells", "monocytes", "eosinophils"],
        "notes":   "Paired with Siglec-9 in suppressive signalling on NK/myeloid cells. "
                   "CD24-Siglec-10 axis is an emerging immune checkpoint in cancer.",
        "doi": [
            "https://doi.org/10.1073/pnas.2026102118",
            "https://doi.org/10.3390/cells10051260",
        ],
    },

    "Siglec-15": {
        "abbreviation":    ["Siglec-15", "CD33L3"],
        "name":    "Siglec-15",
        "gene":    "SIGLEC15",
        "family":  "Siglec",
        "uniprot": "Q6ZMC8",
        "primary": ["Neu5Ac(a2-6)GalNAc", "Sialyl-Tn"],
        "secondary":["Neu5Ac(a2-3)Gal"],
        "negative": [],
        "tissues": ["osteoclasts", "macrophages", "tumor-associated macrophages"],
        "notes":   "Promotes osteoclast fusion and bone resorption; expressed on osteoclast "
                   "precursors. Emerging immune checkpoint in cancer — structurally related "
                   "to PD-L1 in some tumour contexts. Target for antibody-based therapy.",
        "doi": [
            "https://doi.org/10.1038/nm.2834",
            "https://doi.org/10.1016/j.celrep.2019.10.113",
        ],
    },

    "E-Selectin": {
        "abbreviation":    ["E-Selectin", "SELE", "CD62E", "ELAM-1"],
        "name":    "E-Selectin",
        "gene":    "SELE",
        "family":  "Selectin",
        "uniprot": "P16581",
        "primary": ["Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc"],
        "secondary":["Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc",
                     "6-sulfo sialyl Lewis X"],
        "negative": [],
        "tissues": ["activated endothelium", "inflamed venules"],
        "notes":   "Inducible on endothelium by IL-1beta/TNF-alpha. Mediates slow rolling of "
                   "leukocytes. Sialyl Lewis X is the canonical ligand on PSGL-1.",
        "doi": [
            "https://doi.org/10.1182/blood-2011-07-363838",
            "https://doi.org/10.1021/bi0507130",
        ],
    },

    "L-Selectin": {
        "abbreviation":    ["L-Selectin", "SELL", "CD62L", "LAM-1"],
        "name":    "L-Selectin",
        "gene":    "SELL",
        "family":  "Selectin",
        "uniprot": "P14151",
        "primary": ["6-sulfo sialyl Lewis X", "6-sulfo-Gal(b1-4)GlcNAc"],
        "secondary":["Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc"],
        "negative": [],
        "tissues": ["leukocytes", "lymphocytes", "neutrophils"],
        "notes":   "Constitutive on naive leukocytes; mediates homing to lymph nodes via "
                   "PNAd (6-sulfo LacNAc) on HEV. Shed by ADAM17 upon activation.",
        "doi": [
            "https://doi.org/10.1182/blood-2011-07-363838",
            "https://doi.org/10.1021/bi0507130",
        ],
    },

    "P-Selectin": {
        "abbreviation":    ["P-Selectin", "SELP", "CD62P", "GMP-140"],
        "name":    "P-Selectin",
        "gene":    "SELP",
        "family":  "Selectin",
        "uniprot": "P16109",
        "primary": ["Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc",
                    "Core 2 O-glycan sLe^x on PSGL-1"],
        "secondary":["6-sulfo sialyl Lewis X",
                     "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc"],
        "negative": [],
        "tissues": ["platelets (Weibel-Palade bodies)", "activated endothelium"],
        "notes":   "Rapidly translocated from alpha-granules to platelet surface upon activation. "
                   "Canonical ligand is sulfo-tyrosine + Core 2 sLe^x on PSGL-1.",
        "doi": [
            "https://doi.org/10.1182/blood-2011-07-363838",
            "https://doi.org/10.1006/bbrc.2000.3784",
        ],
    },

    "DC-SIGN": {
        "abbreviation":    ["DC-SIGN", "CD209", "CLEC4L"],
        "name":    "DC-SIGN (Dendritic Cell-Specific ICAM-3-Grabbing Non-integrin)",
        "gene":    "CD209",
        "family":  "C-type lectin (CLR)",
        "uniprot": "Q9NNX6",
        "primary": ["Man(a1-2)Man", "Man(a1-3)[Man(a1-6)]Man",
                    "Fuc(a1-2)Gal", "Lewis B", "Lewis Y"],
        "secondary":["GlcNAc(b1-2)Man", "Man(a1-6)Man"],
        "negative": ["sialylated glycans"],
        "tissues": ["dendritic cells", "macrophages (liver, lung, mucosa)"],
        "notes":   "Ca2+-dependent; recognises both high-mannose N-glycans and fucosylated "
                   "blood group antigens. Key pathogen-recognition receptor (HIV, Mycobacteria, "
                   "SARS-CoV-2). Lewis antigens on H. pylori hijack DC-SIGN.",
        "doi": [
            "https://doi.org/10.1016/j.febslet.2006.10.009",
            "https://doi.org/10.1093/glycob/cwq086",
        ],
    },

    "MBL": {
        "abbreviation":    ["MBL", "MBL2", "Mannose-binding lectin"],
        "name":    "Mannose-binding lectin 2",
        "gene":    "MBL2",
        "family":  "C-type lectin (Collectin)",
        "uniprot": "P11226",
        "primary": ["Man(a1-2)Man", "GlcNAc(b1-2)Man", "Fuc(a1-2)Gal"],
        "secondary":["GlcNAc", "Man"],
        "negative": ["Gal", "Neu5Ac"],
        "tissues": ["liver (secreted into serum)", "serum"],
        "notes":   "Serum complement lectin pathway initiator. Trimeric bouquet of CRDs; "
                   "recognises repeating mannose/GlcNAc on pathogen surfaces. "
                   "MBL2 deficiency leads to recurrent infections.",
        "doi": [
            "https://doi.org/10.1093/glycob/cwq086",
            "https://doi.org/10.1093/glycob/21.4.512",
        ],
    },

    "Langerin": {
        "abbreviation":    ["Langerin", "CD207", "CLEC4K"],
        "name":    "Langerin / CD207",
        "gene":    "CLEC4K",
        "family":  "C-type lectin (CLR)",
        "uniprot": "Q9UJ71",
        "primary": ["Man(a1-2)Man", "Fuc(a1-2)Gal", "GlcNAc"],
        "secondary":["Sulfated glycans", "6-sulfo-GlcNAc"],
        "negative": [],
        "tissues": ["Langerhans cells (skin, cornea, mucosa)"],
        "notes":   "Forms Birbeck granules; drives CLR-mediated antigen uptake. "
                   "Selectively restricts HIV infection in Langerhans cells.",
        "doi": [
            "https://doi.org/10.1002/cbic.202000238",
        ],
    },

    "MGL": {
        "abbreviation":    ["MGL", "CD301", "CLEC10A", "ASGR-like"],
        "name":    "Macrophage Galactose-type Lectin (MGL)",
        "gene":    "CLEC10A",
        "family":  "C-type lectin (CLR)",
        "uniprot": "Q8IYK4",
        "primary": ["GalNAc", "Tn antigen (GalNAc-alpha-Ser/Thr)", "Sialyl-Tn"],
        "secondary":["Gal(b1-3)GalNAc (Core 1)", "GalNAc(b1-4)GlcNAc"],
        "negative": ["sialylated Core 1 / Core 2"],
        "tissues": ["dendritic cells", "macrophages (skin, gut, tumour)"],
        "notes":   "Uniquely recognises terminal GalNAc / Tn antigen on tumour mucins "
                   "(MUC1, MUC16). Internalises antigen for cross-presentation. "
                   "Target of cancer immunotherapy.",
        "doi": [
            "https://doi.org/10.1002/cbic.202000238",
            "https://doi.org/10.1038/nrc1251",
        ],
    },

    "ASGR1": {
        "abbreviation":    ["ASGR1", "ASGPR", "Ashwell receptor"],
        "name":    "Asialoglycoprotein receptor 1",
        "gene":    "ASGR1",
        "family":  "C-type lectin (Hepatic lectin)",
        "uniprot": "P07306",
        "primary": ["Gal", "GalNAc"],
        "secondary":["Gal(b1-4)GlcNAc", "GalNAc(b1-4)GlcNAc"],
        "negative": ["Neu5Ac(a2-3/6)Gal"],
        "tissues": ["hepatocytes (liver)", "sinusoidal endothelium"],
        "notes":   "Mediates clearance of desialylated (asialylated) serum glycoproteins. "
                   "Exploited for targeted drug delivery to liver (GalNAc-siRNA conjugates).",
        "doi": [
            "https://doi.org/10.1093/glycob/cwq086",
        ],
    },

    "SP-D": {
        "abbreviation":    ["SP-D", "SFTPD", "Surfactant Protein D"],
        "name":    "Surfactant Protein D",
        "gene":    "SFTPD",
        "family":  "Collectin (C-type lectin)",
        "uniprot": "P35247",
        "primary": ["Man(a1-2)Man", "GlcNAc(b1-2)Man", "Fuc(a1-2)Gal"],
        "secondary":["Man", "GlcNAc"],
        "negative": ["Galactose (minimal binding)", "Sialic acid"],
        "tissues": ["lung — alveolar epithelium (type II pneumocytes)",
                    "alveolar macrophages"],
        "notes":   "Pulmonary collectin; innate immune defense in lung. Forms cruciform "
                   "dodecamers that agglutinate pathogens via carbohydrate recognition. "
                   "Binds influenza hemagglutinin and neutralises many respiratory pathogens.",
        "doi": [
            "https://doi.org/10.1146/annurev.immunol.18.1.541",
            "https://doi.org/10.1016/j.jaci.2008.02.042",
        ],
    },

    "SP-A": {
        "abbreviation":    ["SP-A", "SFTPA1", "Surfactant Protein A"],
        "name":    "Surfactant Protein A",
        "gene":    "SFTPA1",
        "family":  "Collectin (C-type lectin)",
        "uniprot": "Q8IWL2",
        "primary": ["Man", "Fuc(a1-2)Gal", "GalNAc"],
        "secondary":["GlcNAc", "Man(a1-2)Man"],
        "negative": [],
        "tissues": ["lung — alveolar type II cells", "bronchiolar Clara cells"],
        "notes":   "Most abundant pulmonary collectin (surfactant protein). "
                   "Regulates surfactant metabolism and provides innate defence. "
                   "Binds phospholipid head-groups and carbohydrate moieties simultaneously.",
        "doi": [
            "https://doi.org/10.1152/ajplung.00484.2007",
            "https://doi.org/10.1146/annurev.immunol.18.1.541",
        ],
    },

    "Intelectin-1": {
        "abbreviation":    ["Intelectin-1", "ITLN1", "Omentin", "HL-1"],
        "name":    "Intelectin-1 (Omentin-1)",
        "gene":    "ITLN1",
        "family":  "Intelectin (X-type lectin)",
        "uniprot": "Q8WWA0",
        "primary": ["Galactofuranose (Galf)",
                    "D-Glycero-D-talo-oct-2-ulosonic acid (KO)"],
        "secondary":["Ribose", "D-Fructose"],
        "negative": ["Galactopyranose (mammalian Gal — not bound)"],
        "tissues": ["gut — goblet cells of small intestine / colon",
                    "heart — cardiac myocytes (secreted as adipokine from omentum)"],
        "notes":   "Distinguishes microbial galactofuranose (absent in mammals) from mammalian "
                   "galactopyranose — a mechanism for selective pathogen recognition. "
                   "Also known as omentin-1 with insulin-sensitising metabolic functions.",
        "doi": [
            "https://doi.org/10.1016/j.str.2011.09.017",
            "https://doi.org/10.1093/glycob/cwr094",
        ],
    },

    "M6PR": {
        "abbreviation":    ["M6PR", "CI-MPR", "CD222", "IGF2R", "Cation-independent MPR"],
        "name":    "Mannose-6-Phosphate Receptor (CI-MPR)",
        "gene":    "IGF2R",
        "family":  "P-type lectin",
        "uniprot": "P11717",
        "primary": ["Man-6-phosphate (Man-6-P)", "Phosphomannosyl residues"],
        "secondary":["Man(a1-2)Man-6-P"],
        "negative": ["Unphosphorylated mannose"],
        "tissues": ["ubiquitous — all nucleated cells (Golgi / endosome / plasma membrane)",
                    "hepatocytes", "fibroblasts", "macrophages"],
        "notes":   "Fundamental to lysosomal enzyme trafficking — directs newly synthesised "
                   "lysosomal hydrolases from the Golgi to lysosomes via Man-6-P recognition. "
                   "Also binds IGF-II and TGF-beta. Defects cause I-cell disease.",
        "doi": [
            "https://doi.org/10.1146/annurev.biochem.66.1.303",
            "https://doi.org/10.1038/nrm2444",
        ],
    },
}
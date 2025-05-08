def chrom_insulation(clr, chrom, resolution = 1000, n = 2, size = 4):
    if not "weight" in clr.pixels().columns:
        from cooler import balance_cooler
        balance_cooler(clr, cis_only = True, store = True)
    from cooltools import insulation
    ref_arms = chrom_arms(filt = False)
    ref_arms = ref_arms[ref_arms["chrom"].isin(clr.chromnames)]
    ref_arms["name"] = range(len(ref_arms))
    insulation_tab = insulation(clr, [10*resolution], ref_arms, nproc = n)
    insulation_strength = insulation_tab[f"boundary_strength_{10*resolution}"].tolist()
    insulation_strength_sorted = np.sort(np.array(insulation_strength)[~np.isnan(np.array(insulation_strength))])
    total_size = sum([chrom_sub["end"].max() for chr, chrom_sub in ref_arms.groupby("chrom")])
    average_size = total_size//size 
    # 
    insulate_cutoff = insulation_strength_sorted[-size:].min()
    insulated_filtered = insulation_tab[insulation_tab[f"boundary_strength_{10*resolution}"] >= insulate_cutoff]
    return insulation_filtered
    
def domain_bed_parse(domain_group, sig_contact, resolution = 1000):
    """ 
    Return: 
    domain in dataframe of bed format
    """
    pass_col = [col for col in sig_contact.columns if "_pass" in col]
    sig_contact_filtered = sig_contact[sig_contact[pass_col].mean(axis = 1) == 1].copy()
    domain_list = list()
    for a in domain_group:
        sign_domain = sig_contact_filtered[(sig_contact_filtered.bin1.isin(domain_group[a])) | (sig_contact_filtered.bin2.isin(domain_group[a]))].copy()
        sign_domain["domain"] = a
        domain_list.append(sign_domain)
    domain_df = pd.concat(domain_list, axis = 0)
    domain_range = list()
    for (chrom, d), domain_sub in domain_df.groupby(["chrom", "domain"]):
        domain_range.append((chrom, domain_sub["bin1"].min()*resolution, domain_sub["bin2"].max()*resolution))
    domain_bed = pd.DataFrame(domain_range, columns = ["chrom", "start", "end"])
    return domain_bed

def connected_domain(sig_contact):
    """
    Return: 
    dict: each domain store as a key and its bin number
    """
    sig_contact.sort_values(by = ["bin1", "bin2"], inplace = True, ignore_index = True)
    il = sig_contact.iloc[0]["bin1"]; ir = sig_contact.iloc[0]["bin2"]
    a = 0
    domain_group = dict()
    domain_group[a] = [sig_contact.iloc[0].bin1, sig_contact.iloc[0].bin2]
    for row in sig_contact.iloc[1:].itertuples():
        rl = row.bin1
        rr = row.bin2
        if rl > ir:
            a += 1
            domain_group[a] = [rl, rr]
            ir = rr
        if ir >= rl:
            domain_group[a].append(rl); domain_group[a].append(rr)
            ir = max(rr, ir)
    return domain_group


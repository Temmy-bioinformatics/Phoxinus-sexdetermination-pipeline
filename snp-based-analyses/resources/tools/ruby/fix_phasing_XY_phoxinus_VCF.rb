# m_matschiner Thu Mar 26 19:26:33 CET 2020

# Read the vcf file.
invcf_file_name = ARGV[0]
invcf_file = File.open(invcf_file_name)
invcf_lines = invcf_file.readlines

# Define males and females.
###Define males as males for XY and females as males for the ZW species, it will fix according to heterogametic sex

####In this case, females are defined as males because it is a ZW system
# List of females (2 denotes females)
females = [
    "Agger_ZFMK-TIS-65414", "Agger_ZFMK-TIS-65499", "Agger_ZFMK-TIS-66094",
    "Agger_ZFMK-TIS-66100", "Agger_ZFMK-TIS-66114", "Agger_ZFMK-TIS-66357",
    "Agger_ZFMK-TIS-66363", "Anger_ZFMK-TIS-80124", "Anger_ZFMK-TIS-80125",
    "Naafbach_ZFMK-TIS-66031", "Naafbach_ZFMK-TIS-66058", "Suelz_ZFMK-TIS-66172",
    "Suelz_ZFMK-TIS-66201"
]

# List of males (1 denotes males)
males = [
    "Agger_ZFMK-TIS-65468", "Agger_ZFMK-TIS-65486", "Agger_ZFMK-TIS-65487",
    "Agger_ZFMK-TIS-66086", "Naafbach_ZFMK-TIS-66030", "Naafbach_ZFMK-TIS-66046",
    "Naafbach_ZFMK-TIS-66047", "Agger_ZFMK-TIS-66336", "Anger_ZFMK-TIS-80141",
    "Anger_ZFMK-TIS-80142", "Mouzon_ZFMK-TIS-80148", "Mouzon_ZFMK-TIS-80150",
    "Suelz_ZFMK-TIS-66168"
]
##rather heterozygotic ZW or XY sex

out_vcf_string = ""
samples = []
seqs1 = []
seqs2 = []
invcf_lines.each do |l|
	if l[0..1] == "##"
		out_vcf_string << l
	elsif l[0] == "#"
		out_vcf_string << l
		samples = l.split[9..-1]
		samples.size.times do
			seqs1 << []
			seqs2 << []
		end
	else
		# Get the reference and alt allele and sample genotypes.
		line_ary = l.split
		alleles = Marshal.load(Marshal.dump(line_ary[3]))
		line_ary[4].split(",").each do |a|
			alleles << a
		end

		gts = line_ary[9..-1]
		samples.size.times do |x|
			gt = gts[x]
			if males.include?(samples[x])
				female = females[males.index(samples[x])]
				female_gt = gts[samples.index(female)]
				if female_gt.split("|")[0] == female_gt.split("|")[1]
					gt_dist = 0
					gt_dist += 1 if gt.split("|")[0] != female_gt.split("|")[0]
					gt_rev_dist = 0
					gt_rev_dist += 1 if gt.split("|")[1] != female_gt.split("|")[0]
					if gt_rev_dist < gt_dist
						gts[x] = gt.reverse
					end

				end
			end
		end

		# Prepare the phased sequences.
		sample_i = 0
		gts.each do |gt|
			this_gt1 = gt.split("|")[0]
			this_gt2 = gt.split("|")[1]
			allele1 = alleles[this_gt1.to_i]
			allele2 = alleles[this_gt2.to_i]
			allele1 = "N" if allele1 == "*"
			allele2 = "N" if allele2 == "*"
			seqs1[sample_i] << allele1
			seqs2[sample_i] << allele2
			sample_i += 1
		end

		# Add to the vcf output.
		out_vcf_string << line_ary[0..8].join("\t")
		out_vcf_string << "\t"
		out_vcf_string << gts.join("\t")
		out_vcf_string << "\n"

	end
end

# Write the phylip output file.
out_phylip_string = "#{samples.size*2} #{seqs1[0].size}\n"
samples.size.times do |x|
	out_phylip_string << "#{samples[x]}_1  #{seqs1[x].join("")}\n"
	out_phylip_string << "#{samples[x]}_2  #{seqs2[x].join("")}\n"
end
out_phylip_name = invcf_file_name.sub(".vcf",".phy")
out_phylip = File.open(out_phylip_name,"w")
out_phylip.write(out_phylip_string)

# Write the vcf output file.
out_vcf_name = invcf_file_name.sub(".vcf","_fixed.vcf")
out_vcf = File.open(out_vcf_name,"w")
out_vcf.write(out_vcf_string)

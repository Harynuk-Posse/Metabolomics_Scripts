Scripts used for classification of compounds 

Non-derivatized  ‘***************************************************************************
function straightchain_aldehyde() as Boolean
'high abundance m/z41, 55 
Dim Carbon_number
Dim Expected_MW
Dim Em 
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed1 = (Rank(1) = 41 or Rank(1) = 44 or rank(1)=57) 
bed2 = abundance(41)>600 and abundance(42)<500 and abundance(44)>200 and abundance(45)>50
bed3 = abundance(55)>150 and abundance(57)>400
bed4 = abundance(74)<50
bed5 = (abundance(57)/abundance(55)>0.95) 
bed6 = abundance(82)>abundance(83)

ionseriescounter1 = 0 
for k = 47 to 52 step 1 
	if abundance(k)>60 then
		ionseriescounter1 = ionseriescounter1+1
	end if 
next
ionseriescounter2 = 0 
for m = 74 to 79 step 1 
	if abundance(m)>50 then
		ionseriescounter2 = ionseriescounter2+1
	end if 
next
ionseriescounter3 = 0 
for n = 87 to 94 step 1 
	if abundance(n)>60 then
		ionseriescounter3 = ionseriescounter3+1
	end if 
next 

Carbon_number = 4
Do while Carbon_number<=30
Expected_MW = 14*Carbon_number +16
	m18 = Expected_MW -18
	m28 = Expected_MW -28
	m43 = Expected_MW -43
	m44 = Expected_MW -44
	noise_counter = 0
	for noisecheck = Expected_MW+2 to Em step 1 
		If Intensity(noisecheck) <= (ave_intensity + 6*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter = noise_counter + 1
		end if	
next	
If noise_counter < 5 and abundance(m28)>0 and abundance(m43)>0 and abundance(m44)>10 and ionseriescounter1<=1 and ionseriescounter2<=1 and ionseriescounter3<=1 and bed1 and bed2 and bed3 and bed4 and bed5 and bed6 then straightchain_aldehyde =true  
Carbon_number = Carbon_number+1
Loop
end function
‘***************************************************************************
Function Ketones()as Boolean
'base peak 43
'high abundance peak ion series m/z 71,85,99 
Dim Carbon_number
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity
Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5

bed1 = rank(1) = 43 or rank(1) = 58
bed2 = abundance(43)>700 and abundance(71)>50 and abundance(57)<200 and abundance(59)>50 and abundance(58)>200

fragment_counter = 0
for n = 5 to 30  step 1
	fragment_test =n*14+1
	if abundance(fragment_test)>10 then
		fragment_counter= fragment_counter +1	
	end if
next	

Carbon_number = 5
Do while Carbon_number<=30 
Expected_MW = 14*Carbon_number+16
M43 = Expected_MW - 43
noise_counter1 = 0
noise_counter2 = 0 
	for noisecheck = Expected_MW+2 to Em step 1 
		If Intensity(noisecheck) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and abundance(M43)>0 and fragment_counter >= (Carbon_number-5) and noise_counter1 <=5 and noise_counter2<5 then Ketones = true
Carbon_number = Carbon_number+1
Loop	
end function 
‘***************************************************************************
function alcohols2() as Boolean
'2-alcohol 
‘base peak m/z45
'ion series 55,69,83,97,111,125
'aldehyde 'ion series 68,82,96,110,124,138 

Dim Carbon_number
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity
Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5

bed1 = Rank(1)=45
bed2 = abundance(41)>100 and abundance(43)>100 
bed3 = abundance(55)>100 and abundance(56)>30 and abundance(57)>10

fragment_counter = 0
for n = 3 to 30  step 1
	fragment_alcohol =n*14-1
	if abundance(fragment_alcohol)>0 then
		fragment_counter= fragment_counter +1	
	end if
next	

Carbon_number = 4 
Do while Carbon_number<=30
	Expected_MW = 14*Carbon_number+18		
	m18 = Expected_MW -18
	noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
	noise_counter3 = 0
	for noisecheck2 = 100 to Em step 1 
		If abundance(noisecheck2)> 100 then
			noise_counter3 = noise_counter3 + 1
		end if
	next
If bed1 and bed2 and bed3 and abundance(m18)>5 and fragment_counter >= (Carbon_number-5) and noise_counter1<=5 and noise_counter2<5 and noise_counter3<1 then Alcohols2 = true
Carbon_number = Carbon_number+1
Loop	
End function
‘***************************************************************************
function 1-alcohols() as Boolean
'alcohol  'ion series 55,69,83,97,111,125
'aldehyde 'ion series 68,82,96,110,124,138 

Dim Carbon_number
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity
Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed1 = Rank(1)=55 or Rank(1)=43 or Rank(1)=41
bed2 = abundance(55)>700 and abundance(56)>400 and abundance(57)>200
bed3 = abundance(69)>400 and abundance(83)>200

fragment_counter = 0
for n = 3 to 30  step 1
	fragment_alcohol =n*14-1
	fragment_aldehyde =n*14-2
	fragment_next = n*14+1
	if abundance(fragment_alcohol)>abundance(fragment_aldehyde) and abundance(fragment_alcohol)>abundance(fragment_next) and abundance(fragment_alcohol)>10 then
		fragment_counter= fragment_counter +1	
	end if
next	
Carbon_number = 4 
Do while Carbon_number<=30 
	Expected_MW = 14*Carbon_number+18		
	m18 = Expected_MW -18
	m46 = Expected_MW-46
	noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
	noise_counter3 = 0
	for noisecheck2 = 100 to Em step 1 
		If abundance(noisecheck2)> 100 then
			noise_counter3 = noise_counter3 + 1
		end if
	next
If bed1 and bed2 and bed3 and abundance(m18)>0 and abundance(m46)>10 and fragment_counter >= (Carbon_number-5) and noise_counter1<=5 and noise_counter2<5 and noise_counter3<1 then 1-alcohols = true
Carbon_number = Carbon_number+1
Loop
End function
‘***************************************************************************
function linearSC_FFA() as Boolean
'[HOOC(CH2)n]+ from m/z = 115 to 255
'Checks for ion series of general formula [CH3OCO(CH2)n]+ 
'base peak at m/z=60
'high abundance of m/z73
'[M-17] loss of OH-
'[M-29] & [M-43]

Dim Carbon_number
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity
 
Em = EndMass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed1 = (Rank(1)=60 and Rank(2)=73)
bed2 = abundance(60)>500 and abundance(61)>10 and abundance(73)>300 and abundance(55)<500
bed3 = abundance(41)>100 and abundance(42)>50 and abundance(43)>50 and abundance(45)>100
Carbon_number = 4
Do while Carbon_number < 15
	Expected_MW = 14*(Carbon_number-1)+46
	m17 = Expected_MW -17		
m29 = Expected_MW -29
	m43 = Expected_MW -43
	bed4 = abundance(m17)>0 or abundance(m29)>0 or abundance(m43)>0
noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
	if bed1 and bed2 and bed3 and bed4 and noise_counter1<=5 and noise_counter2<5 then linearSC_FFA =true			
	Carbon_number = Carbon_number + 1
	Loop
end function
‘***************************************************************************
function linearLC_FFA()as Boolean
Dim Carbon_number
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = EndMass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed1 = (Rank(1)=60 or Rank(1)=73)
bed2 = Abundance(73)>700 and abundance(60)>700
bed3 = abundance(55)>400 and abundance(57)>200
bed4 = abundance(83)>50 and abundance(85)>30 and abundance(87)>50 and abundance(129)>100 
bed5 = abundance(41)>500 and abundance(43)>500
bed6 = abundance(97)>10 and abundance(61)>100 and abundance(61)<400 and abundance(69)>100

'ion series [HOOC(CH2)n]+ from m/z = 115 to 255
fragment_counter = 0
for n = 8 to 30  step 1
	fragment_test =n*14+3
	if abundance(fragment_test)>5 then
		fragment_counter= fragment_counter +1	
	end if
next	
Carbon_number = 8
Do while Carbon_number < 30
	Expected_MW = 14*(Carbon_number-1)+46
	m29 = Expected_MW -29
	m43 = Expected_MW -43
noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
if abundance(Expected_MW)>5 and abundance(m29)>0 and abundance(m43)>0 and bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and noise_counter1<=5 and noise_counter2<5 and fragment_counter>=Carbon_number-9 then linearLC_FFA =true		
Carbon_number = Carbon_number + 1
Loop
end function
‘***************************************************************************
'm/z 43 base peak 
'm/z 60 rank(2)
'high abundance(102)
'ion series 45,59,73,87,101,115,129,143,157,171  
function isopropylester() as Boolean
Dim Carbon_number
Dim Expected_MW
Dim Em 
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed3 = abundance(102)>300 and abundance(97)>50
bed4 = Rank(1)=43 
bed5 = (Rank(2)=60 or Rank(3)=60)

fragment_counter1 = 0
for n = 3 to 30  step 1
	fragment_test =n*14+3
	if abundance(fragment_test)>15 then
		fragment_counter1= fragment_counter1+1	
	end if
next	

Carbon_number = 4
Do while Carbon_number < 30
	Expected_MW = 14*Carbon_number+74
	m40 = Expected_MW -40
	m41 = Expected_MW -41
	m42 = Expected_MW -42
	m59 = Expected_MW -59
	bed1 = abundance(m40)>1 and abundance(m41)>30 and abundance(m42)>50
bed2 = abundance(m59)>30
noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and fragment_counter1>=(Carbon_number-6) and noise_counter1<=5 and noise_counter2<5 then isopropylester = true
Carbon_number = Carbon_number + 1
Loop	
End function
‘***************************************************************************
'm/z 88 base peak instead of m/z 74
‘high abundance m/z 101
‘high abundance of [M-45] & [M-43] (loss of the ethoxide)
‘[M-29] loss of the ethyl group 
'ion series 101,115,129,143,157,171,185,199 
‘m/z 41,55,60,73 
function FA_EthylEsters() as Boolean
Dim Carbon_number
Dim Expected_MW
Dim Em 
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed1 = Rank(1)=88 
bed2 = abundance(41)>200 and abundance(43)>200
bed3 = abundance(55)>100 and abundance(57)>100 and abundance(60)>50 and abundance(61)>50
bed4 = abundance(73)>100 and abundance(101)>200 

fragment_counter = 0
for n = 3 to 30  step 1
	fragment_test =n*14+3
	if abundance(fragment_test)>3 then
		fragment_counter= fragment_counter+1	
	end if
next	

Carbon_number = 4
Do while Carbon_number < 30
	Expected_MW = 14*Carbon_number+60
m27 = Expected_MW -27
m29 = Expected_MW -29
	m43 = Expected_MW -43
	m45 = Expected_MW -45
bed5 = abundance(m27)>0 or abundance(m29)>0  
bed6 = abundance(Expected_MW)>0 and abundance(m43)>3 and abundance(m45)>3
noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and fragment_counter>=(Carbon_number-6) and noise_counter1<=5 and noise_counter2<5 then FA_EthylEsters = true
Carbon_number = Carbon_number + 1
Loop	
End function
‘***************************************************************************
Function LinearSaturatedFAME_v2() as boolean

Dim Carbon_number
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5

bed1 = Rank(1)=74 
bed2 = Rank(2)=87 or Rank(3)=87
bed3 = (Abundance(55)>150 And Abundance(55)<400) 
bed4 = abundance(57)<350
bed5 = abundance(59)>50

'Checks for ion series of general formula [CH3OCO(CH2)n]+ 

fragment_counter = 0
for n = 2 to 30  step 1
	fragment_test =n*14+59
	if abundance(fragment_test)>5 then
		fragment_counter= fragment_counter +1	
	else
		LinearSaturatedFAME_v2 = false		
	end if
next

Carbon_number = 5 
Do while Carbon_number < 30	 
Expected_MW = 14*(Carbon_number-1) + 60
M29 = Expected_MW -29
M31 = Expected_MW -31
M43 = Expected_MW -43
bed6 = Abundance(Expected_MW)>5 
bed7 = Abundance(M29)>0 or Abundance(M31)>0 or Abundance(M43)>0 
noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and fragment_counter >= (Carbon_number-10) and noise_counter1<=5 and noise_counter2<5 then LinearSaturatedFAME_v2 = true
Carbon_number = Carbon_number + 1
Loop	
end function

'***************************************************************************
'Linear Monoenoic Fatty Acids Methyl Esters
'1. base peak: rank(1) = 55
'2. m/z69 in high abundance, compare to saturated FAMEs 
'm/z69>300
'3. molecular ion at reasonable abundance
'4. m/z74in lower abundance, compare to saturated FAMEs 
'm/z74<600
'5. m/z83 in higher abundance, compare to saturated FAMEs 
'm/z83>200
'6. long homologous series (general formula [CnH2n 1]+)  of ions 14 amu apart at m/z = 83, 97, 111, 125, 139, 153, 167
'7. high abundance: [M-32] & [M-74]& [M-116]
'8. Presence of [M-60] & [M-61] and [M-49] & [M-50]
'***************************************************************************
function Linear_monoenoicFAME_v2() as boolean

Dim Carbon_number
Dim Expected_MW
Dim Em 
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5

bed1 = Rank(1)=55
bed2 = Abundance(69)>300    
bed3 = abundance(74)<600
bed4 = abundance(83)>150
bed5 = abundance(74)/abundance(69)>0.8 and abundance(74)/abundance(69)<1.6

'Checks for ion series of general formula [CnH2n-1]+) 
'Prominent features and must have fragments

fragment_counter = 0
for n = 6 to 30  step 1
	fragment_test =n*14- 1
	if abundance(fragment_test)>5 then
		fragment_counter= fragment_counter +1	
	end if
next

Carbon_number = 5 
Do while Carbon_number < 30
Expected_MW = 14*(Carbon_number-1) + 58
M32 = Expected_MW -32
M74 = Expected_MW -74
M116 = Expected_MW -116
M49 = Expected_MW -49 
M50 = Expected_MW -50
M60 = Expected_MW -60 
M61 = Expected_MW -61

bed6 = Abundance(M32)>10 and Abundance(M74)>10 and Abundance(M116)>10
bed7 =  Abundance(M49)>0 and Abundance(M50)>0 
bed8 = Abundance(M60)>0 and Abundance(M61)>0 
noise_counter1 = 0
noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and bed8 and fragment_counter >= (Carbon_number-8) and noise_counter1<=5 and noise_counter2<5 then Linear_monoenoicFAME_v2= true
Carbon_number = Carbon_number + 1
Loop	
End function
'***************************************************************************
function Linear_dienoicFAME_v2() as boolean

Dim Carbon_number
Dim Expected_MW
Dim Em 
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity
Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5

bed1 = Rank(1)=67 
bed2 = Abundance(74)<250
bed3 = Abundance(81)>500 and abundance(81)>abundance(79)     
bed4 = abundance(68)<500 and abundance(82)<500 and abundance(73)<200 and abundance(75)<200
bed5 = (abundance(95)>abundance(96)) and (abundance(109)>abundance(110))
'Checks for Ion Series CnH2n-3 
'm/z 67,81,95,109,123,
fragment_counter = 0
for n = 5 to 30 step 1
	fragment_test1 = n*14-3
	if abundance(fragment_test1)>0 then
		fragment_counter =  fragment_counter +1
	end if
next
'prominent features and must have fragments
'base peak 67 or 81
Carbon_number = 5
Do while Carbon_number < 30
Expected_MW = 14*(Carbon_number-1) + 56
M31 = Expected_MW -31
M32 = Expected_MW -32
bed6 = abundance(M31)>0 and Abundance(M32)>0 
noise_counter1 = 0
noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and fragment_counter >= (Carbon_number-6) and abundance(Expected_MW)>0 and noise_counter1<=5 and noise_counter2<5 then Linear_dienoicFAME_v2= true
Carbon_number = Carbon_number + 1
Loop	
end function
'***************************************************************************
function Linear_trienoicFAME_v2() as Boolean
Dim Carbon_number
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5 
bed1 = Rank(1)= 79
bed2 = Abundance(74)<250 and abundance(59)>100
ratio91=abundance(91)/abundance(93)
bed3 = ratio91<1.3
bed4 = (Rank(2)=41 or Rank(2)=55 or Rank(2) =67) 
bed5 = abundance(55)>400 and abundance(67)>500
bed6 = abundance(95)>100 and abundance(108)>30 and abundance(121)>50 and abundance(135)>30
'ion series m/z 65,79,93,107,121,135,149,163,177
fragment_counter = 0
for n = 5 to 30 step 1
	fragment_test = n*14-5
	if abundance(fragment_test)>0 then
		fragment_counter =  fragment_counter +1
	end if
next
Carbon_number = 5
Do while Carbon_number < 30	 
Expected_MW = 14*(Carbon_number-1) + 54
'M29 = Expected_MW -29
M31 = Expected_MW -31
bed7 = Abundance(M31)>0 
noise_counter1 = 0
noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and fragment_counter >= (Carbon_number-6) and abundance(Expected_MW)>0 and noise_counter1<=5 and noise_counter2<5 then Linear_trienoicFAME_v2= true
Carbon_number = Carbon_number + 1
Loop	
end function

'***************************************************************************
function Linear_multienoicFAME_v2() as Boolean
Dim Carbon_number
Dim Expected_MW
Dim Em 
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed1 = Rank(1)= 79
bed2 = (rank(2)=91 or rank(3)=91)
bed3 = abundance(105)>200 and abundance(105)<400
bed4 = abundance(67)>500 and abundance(67)<900
bed5 = abundance(87)<100 and abundance(55)>200 
bed6 = (abundance(105)/abundance(108))>2.5
bed7 = abundance(74)>50 and abundance(77)>100 and abundance(78)>100 and abundance(80)>100 and abundance(81)>100
'ion series m/z 75,89,103,117,131,145,159,173
fragment_counter = 0
for n = 5 to 30 step 1
	fragment_test = n*14+5
	if abundance(fragment_test)>0 then
		fragment_counter =  fragment_counter +1
	end if
next

Carbon_number = 5
Do while Carbon_number < 30
	Expected_MW = 14*Carbon_number+60
m27 = Expected_MW -27
m29 = Expected_MW -29
	m43 = Expected_MW -43
	m45 = Expected_MW -45
bed5 = abundance(m27)>0 or abundance(m29)>0  
bed6 = abundance(Expected_MW)>0 and abundance(m43)>3 and abundance(m45)>3
noise_counter1 = 0
	noise_counter2 = 0 
	for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and fragment_counter>5 then Linear_multienoicFAME_v2= true	
end function
‘***************************************************************************
TMS derivatized 
 ‘**************************************************************************
function Valine_2TMS_v2() as Boolean
'molecular weight 261
mass = 261
M15 = mass-15
M43 = mass-43
bed1 = (Rank(1)=73 and Rank(2)=144)
bed2 = abundance(M43)>100 and abundance(M15)>0 
bed3 = abundance(73)>500 and abundance(74)>10 and abundance(75)>10
bed4 = abundance(144)>400 and abundance(145)>10 and abundance(147)>10
bed5 = abundance(45)>0 and abundance(45)<400
bed6 = abundance(59)>0 and abundance(59)<250
bed7 = abundance(85)>0 and abundance(100)>10 and abundance(128)>0
bed8 = abundance(218)>30 and abundance(219)>0
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and bed8 then Valine_2TMS_v2 = true
end function
‘***************************************************************************
function proline_2TMS() as boolean
Dim MW
Dim M15
Dim M29
Dim M43

MW=259
M15 = MW-15
M29 = MW-29
M43 = MW -43
bed1 = Rank(1)=73 
bed2 = Rank(2)=142
bed3 = abundance(MW)>0 and abundance(M15)>0 and abundance(M29)>0 and abundance(M43)>30
bed4 = abundance(217)>0 and abundance(218)>0 
bed5 = abundance(84)>5 and abundance(100)>10 and abundance(147)>10 and abundance(170)>0  
bed6 = abundance(45)>100 and abundance(59)>10 and abundance(66)>10 and abundance(75)>10
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 then proline_2TMS = true
end function

‘***************************************************************************
function oxoproline_2TMS() as Boolean
‘MW = 273
bed1 = Rank(1)=73 or Rank(2)= 73
bed2 = Rank(1)=156 or Rank(2)=156
bed3 = abundance(273)>0 and abundance(258)>10 and abundance(230)>10
bed4 = abundance(147)>50 and abundance(133)>5 and abundance(156)>500 
bed5 = abundance(45)>100 and abundance(59)>10 and abundance(84)>3 and abundance(86)>3 and abundance(100)>3
If bed1 and bed2 and bed3 and bed4 and bed5 then oxoproline_2TMS = true
end function

‘***************************************************************************
function tyrosine() as boolean

bed1 = Rank(1)=73 
bed2 = Rank(2)=218 
bed3 = abundance(218)>400 and abundance(219)>50 and abundance(220)>10 and abundance(179)>30
bed4 = abundance(100)>100 and abundance(147)>50  
bed5 = abundance(45)>100 and abundance(59)>10  
If bed1 and bed2 and bed3 and bed4 and bed5 then tyrosine_3TMS = true

bed6 = Rank(1)=73 
bed7 = Rank(2)=179
bed8 = abundance(179)>500 and abundance(180)>100 and abundance(181)>30 and abundance(182)>0
bed9 = abundance(208)>50 and abundance(219)>30  
If bed6 and bed7 and bed8 and bed9 and bed5 then tyrosine_2TMS = true

If tyrosine_3TMS or tyrosine_2TMS then tyrosine = true
end function
'***************************************************************************
Function Lysine_4TMS() as boolean
Dim Em
Em = Endmass()
'TMS-NH=CH-COOTMS m/z 218
bed1 = abundance(434)>0 and abundance(435)>0 'M & [M+1]
bed2 = abundance(329)>0          'M-105
bed3 = abundance(59)>100 and abundance(86)>50 and abundance(100)>100 and abundance(115)>0 and abundance(128)>50 and abundance(147)>0  
bed4 = (Rank(1)=73)
bed5 = (Rank(2)=156 or Rank(2)=174) 
bed6 = (Rank(3)=156 or Rank(3)=174)
bed7 = abundance(186)>0 and abundance(200)>0 and abundance(218)>0 and abundance(230)>10 and abundance(273)>0 and abundance(317)>30  
noise_counter = 0
for noisecheck = 435 to Em step 1 
	If abundance(noisecheck)>10 then
		noise_counter = noise_counter + 1		
	end if 
next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and noise_counter<5 then Lysine_4TMS=true
End function
‘**************************************************************************
‘M+ = 276
Function Asparagine_2TMS() as Boolean
Dim Expected_MW
Dim Em
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity
Expected_MW = 276
Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5

bed1 = Rank(1)=44
bed2 = (Rank(2)=73 or Rank(2)=75 or Rank(2)=159)
bed3 = abundance(73)>200 and abundance(75)>200 
bed4 = abundance(86)>5 and abundance(100>10 and abundance(116)>50       
bed5 = abundance(130)>20 and abundance(147)>20 and abundance(159)>100 and abundance(186)>3 and abundance(244)>3
for noisecheck1 = Expected_MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and noise_counter1<=5 and noise_counter2<5 then Asparagine_2TMS = true
End function
‘***************************************************************************
Function glycine_2TMS 
Dim MW_MonoTMS_glycine 
Dim MW_TriTMS_glycine 
Dim Em
Em = Endmass()
bed1 = (Rank(1)=73 or Rank(2)=73)
bed2 = (Rank(1)=102 or Rank(2)=102)
bed2 = abundance(45)>100 and abundance(59)>100 and abundance(86)>100 and abundance(100)>100 and abundance(133)>50 and abundance(159)>0
bed3 = abundance(174)>400
bed4 = abundance(117)>0 and abundance(147)>200
bed5 = abundance(276)>10 'M-15
bed6 = abundance(248)>50 'M-43

'not significant picks from m/z179 to m/z 245 
significant_counter1 = 0
for sig_ion_check1 = 179 to 245 step 1 
	If abundance(sig_ion_check1)>10 then
		significant_counter1 = significant_counter1 + 1		
	end if
next

significant_counter2 = 0
for sig_ion_check2 = 252 to 275 step 1 
	If abundance(sig_ion_check2)>5 then
		significant_counter2 = significant_counter2 + 1		
	end if
next

noise_counter = 0
for noisecheck = MW_TriTMS_glycine+1 to Em step 1 
	If abundance(noisecheck)>10 then
		noise_counter = noise_counter + 1		
	end if
next

If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and significant_counter1<3 and significant_counter2<2 and noise_counter<5 then TriTMS_glycine=true
End Function
‘***************************************************************************
function MonoenoicFA_TMS() as boolean
Dim Carbon_number
Dim Expected_MW
'ion series m/z 
bed1 = (Rank(1)=75 or rank(2)=75)
bed2 = Abundance(73)>500    
bed3 = abundance(117)>300 and abundance(129)>200 and abundance(145)>50
bed4 = abundance(41)>300 and abundance(55)>300
bed5 = abundance(74)<150 and abundance(91)<100
bed6 = abundance(81)>100 and abundance(96)>100 and abundance(110)>20

'from m/z145 to M-15 
'ion series m/z 157,171,185,199,213,227,241,255,269,283,297 etc
fragment_counter1 = 0
for n = 11 to 30  step 1
	fragment_test =n*14+ 3
	if abundance(fragment_test)>1 then
		fragment_counter1= fragment_counter1 +1	
	end if
next

for Carbon_number = 5 to 30 step 1	 
Expected_MW = 14*(Carbon_number) + 102
M15 = Expected_MW -15 
M16 = Expected_MW -16
bed7 = abundance(Expected_MW)>0 and Abundance(M15)>30 
fragment_counter2 = 0
	for k = 146 to M16
		if abundance(k)>50 then
			fragment_counter2= fragment_counter2 +1	
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and fragment_counter1 >= (Carbon_number-10) and fragment_counter2<3 then MonoenoicFA_TMS = true
next
end function

‘***************************************************************************	
function DienoicFA_TMS() as boolean
Dim Carbon_number
Dim Expected_MW
for Carbon_number = 5 to 30 step 1	 
Expected_MW = 14*(Carbon_number) + 100
M15 = Expected_MW -15 
M16 = Expected_MW -16
bed1 = (abundance(M15)/abundance(Expected_MW))>3
bed2 = (Rank(1)=73 or Rank(1)=75 or Rank(2)=73 or Rank(2)=75)
bed3 = abundance(117)>100 and abundance(129)>100 and abundance(132)<100
bed4 = abundance(41)>200 and abundance(55)>200
bed5 = abundance(67)>200
bed6 = abundance(73)>500 and abundance(75)>500 and abundance(79)<500
bed7 = Abundance(M15)>20 and Abundance(Expected_MW)>0 
bed8 = abundance(150)>3 and abundance(164)>3 
counter1 = 0
	for k=130 to M16 step 1
		if abundance(k)>100 then 
			counter1 = counter1+1 
		end if 
	next 	
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and bed8 and counter1<3 then DienoicFA_TMS = true
next
end function

'***************************************************************************
'multidouble bonds in aliphatic chain
'high abundance of m/z 79

function multienoicFA_TMS() as boolean
Dim Em 
Em = Endmass()
bed1 = (Rank(1)=73 or Rank(1)=75 or rank(1)=79)
bed2 = abundance(79)>500 
bed3 = abundance(108)>10 and abundance(129)>50 and abundance(135)>10
bed4 = abundance(41)>300 and abundance(55)>200 and abundance(67)>400
bed5 = abundance(91)>150 and abundance(93)>100 and abundance(95)>50
bed6 = abundance(73)>500 and abundance(75)>500 
bed7 = abundance(105)>50 and abundance(80)>100 and abundance(81)>30
bed8 = abundance(107)>10 and abundance(108)>10 
counter1 = 0
	for k=160 to Em step 1
		if abundance(k)>100 then 
			counter1 = counter1+1 
		end if 
	next 	
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and bed8 and counter1<3 then multienoicFA_TMS = true
end function
'***************************************************************************
function SatFA_TMS() as boolean
Dim Carbon_number
Dim Expected_MW
Dim Em 
Em = Endmass()
for Carbon_number = 5 to 30 step 1	 
Expected_MW = 14*(Carbon_number) + 104
M16 = Expected_MW -16 	
M15 = Expected_MW -15 	
M14 = Expected_MW -14
M13 = Expected_MW -13
bed1 = abundance(M15)>50 and abundance(M14)>0 and abundance(M13)>0
bed2 = (Rank(1)=73 or Rank(1)=75 or Rank(1)=117)
bed3 = abundance(129)>50 and abundance(130)>0 and abundance(131)>50 and abundance(132)>50
bed4 = abundance(41)>50 and abundance(45)>50 and abundance(55)>50
bed5 = abundance(145)>20
bed6 = abundance(73)>500 and abundance(75)>500 and abundance(117)>200
bed7 = abundance(159)>3 
counter1 = 0
	for n=76 to 116 step 1
		if abundance(n)>90 then 
			counter1 = counter1+1 
		end if 
	next 
counter2 = 0
	for k=146 to M16 step 1
		if abundance(k)>50 then 
			counter2 = counter2+1 
		end if 
	next 	
noise_counter = 0
for noisecheck = Expected_MW+1 to Em step 1 
	If abundance(noisecheck)>15 then
		noise_counter = noise_counter + 1		
	end if 
next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and counter1<3 and counter2<3 and noise_counter<5 then SatFA_TMS = true
next
end function
‘***************************************************************************
Function sugars_TMS()as Boolean
bed1 = rank(1) = 73 
bed2 = abundance(59)>0
bed3 = abundance(147)>0 and abundance(160)>0
bed4 = abundance(103)>50 and abundance(117)>0 and abundance(129)>0 and abundance(133)>0 
bed5 = abundance(204)>0 and abundance(205)>50 and abundance(217)>50 and abundance(229)>0
bed6 = abundance(319)>100 and abundance(320)>0 
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 then sugars_TMS = true
End Function
'***************************************************************************
Function sugar_4TMS() as boolean
Dim Em
Em = Endmass()
bed1 = (Rank(1)=73 and Rank(2)=103)
bed2 = abundance(45)>50 and abundance(59)>30
bed3 = abundance(89)>10 and abundance(103)>300 and abundance(117)>10 and abundance(133)>10 and abundance(147)>50 
bed4 = abundance(160)>10 and abundance(172)>0 and abundance(189)>10 
bed5 = abundance(204)>0 and abundance(205)>0 and abundance(217)>50 and abundance(233)>0
bed6 = abundance(262)>0 and abundance(277)>0 and abundance(307)>10 
bed7 = abundance(74)>0 and abundance(75)>0

noise_counter = 0
for noisecheck = 401 to Em step 1 
	If abundance(noisecheck)>10 then
		noise_counter = noise_counter + 1		
	end if 
next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and noise_counter<5 then sugar_4TMS=true
End function

'***************************************************************************
Function sugar_8TMS() as boolean
Dim Em
Em = Endmass()
bed1 = Rank(1)=73 
bed2 = abundance(45)>30 and abundance(59)>10
bed3 = abundance(81)>5 and abundance(89)>5 and abundance(103)>50  
bed4 = abundance(117)>20 and abundance(129)>50 and abundance(131)>5 and abundance(133)>10 and abundance(147)>100 
bed5 = abundance(157)>10 and abundance(169)>20 and abundance(191)>10 
bed6 = abundance(204)>5 and abundance(205)>0 and abundance(217)>50 and abundance(231)>0 and abundance(243)>10
bed7 = abundance(259)>0 and abundance(271)>10 and abundance(319)>5 and abundance(361)>30
noise_counter = 0
for noisecheck = 401 to Em step 1 
	If abundance(noisecheck)>10 then
		noise_counter = noise_counter + 1		
	end if 
next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and noise_counter<5 then sugar_8TMS=true
End function

'***************************************************************************
Function sugar_5TMS()as Boolean
Dim Em
Em = Endmass()
bed1 = rank(1) = 73 
bed2 = abundance(45)>30 and abundance(59)>10 and abundance(89)>30
bed3 = abundance(103)>50 and abundance(117)>30 and abundance(129)>20 and abundance(131)>5 and abundance(133)>20 
bed4 = abundance(147)>100 and abundance(189)>10 
bed5 = abundance(204)>0 and abundance(205)>10 and abundance(217)>10 and abundance(229)>0
bed6 = abundance(307)>30 
bed7 = abundance(319)>50 and abundance(320)>5 and abundance(320)>1
bed8 = bed6 or bed7
noise_counter = 0
for noisecheck = 465 to Em step 1 
	If abundance(noisecheck)>10 then
		noise_counter = noise_counter + 1		
	end if 
next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed8 and noise_counter<5 then sugar_5TMS=true
End function
‘***************************************************************************
function sterol_TMS() as Boolean
'MW=458
'm/z 41,43,55,57,,73,75
'm/z 81,95,105,119,129,145,159,173,185,199,213,255,275,329,354,368,443,458
'high abundance of m/z 73,129
Dim Carbon_number
Dim Expected_MW
Dim Em 
Dim masscheck1
Dim masscheck2
Dim avg_intensity
Dim sum_intensity
Dim sum_sq
Dim set_intensity

Em = Endmass()
sum_intensity = 0 
noisecounter1 = 0 
for masscheck1 = 500 to Em step 1 'get average noise from m/z 500 to Endnum
	if intensity(masscheck1)>0 then
		sum_intensity = sum_intensity + intensity(masscheck1) 
		noisecounter1 = noisecounter1 + 1	
	end if  
next 
avg_intensity = sum_intensity/noisecounter1
sum_sq = 0
for masscheck2 = 500 to Em step 1
	if intensity(masscheck2)>0 then 
		sq = (intensity(masscheck2)-avg_intensity)^2 
		sum_sq = sum_sq + sq  
	end if
next
stdev_intensity =(sum_sq/(noisecounter1-1))^0.5
bed2 = abundance(41)>200 and abundance(43)>300 and abundance(55)>100 and abundance(57)>50
bed3 = abundance(73)>300 or abundance(75)>300 
bed4 = abundance(91)>100 and abundance(105)>50 and abundance(119)>10 
bed5 = abundance(129)>50 and abundance(145)>10 and abundance(159)>10 and abundance(173)>10 
bed6 = abundance(215)>0 and abundance(233)>0
bed7 = abundance(247)>3 and abundance(255)>10

MW=456 
Do while MW<=500 
noise_counter1 = 0
noise_counter2 = 0 
bed1 = abundance(MW)>0 and abundance(MW-15)>0 and abundance(MW-90)>10 and abundance(105)>5
	for noisecheck1 = MW + 2 to Em step 1 
		If Intensity(noisecheck1) <= (ave_intensity + 4*stdev_intensity) then 
			set_intensity = 0 
		else 
			set_intensity = intensity(noisecheck1) 
		end if 
		If set_intensity/intensity(Rank(1))> 0.005 then
		noise_counter1 = noise_counter1 + 1
		else if abundance(noisecheck1)>30 then
		noise_counter2 = noise_counter2 + 1
		end if
	next
If bed1 and bed2 and bed3 and bed4 and bed5 and bed6 and bed7 and noise_counter1<=5 and noise_counter2<5 then sterol_TMS=true
MW = MW+2
Loop
End function


PRO rm_bkg

;;;;UDPATE HERE;;;;
;useful directories
rawdir = '...' ; direcotry where the calibrated images live
clndir = '...' ; direcotry where the cleaned images should output to
cdedir = '...' ; directory where the code lives

camera = '1';what camera?
ccd = '1';what ccd?
sector = 's0003';what sector?
debug = 1 ;change to 0/1 if you want to write out the fits files
grad1 = 1 ;change to 0/1 if you want to subtract the first gradient
strap = 1 ;change to 0/1 if you want to remove the straps 
grad2 = 1 ;change to 0/1 if you want to subtract any remaining gradient

;;;END UPDATE;;;

bin_box = 64 ;decide on the bins x bins box to find the background, I usally do 64 x 64
bin_itr = 8 ;decide on the size of the iteration to move the box in pixels, I usually do 8
bin_stp = 128 ;decide the size of the pixels used to find the strap level

bin_ex = intarr(bin_box/bin_itr-1) ;make the vector to hold the movement of the box in pixels
jj = 0 ;initialize the index
for ii = bin_itr, bin_box-bin_itr, bin_itr do begin
	bin_ex[jj] = ii 
	jj = jj+1
endfor

;change into the image directory, and get the files to clean
cd, rawdir ; change to raw direcotry
spawn, 'ls *'+sector+'-'+camera+'-'+ccd+'*.fits.gz', files ;get the image list
nfiles = n_elements(files) ;count how many images exist
cd, cdedir; change back to code directory

;loop over each image to clean
for pp = 0l, nfiles-1 do begin
	print, 'Now starting to clean '+files[pp]+' at '+systime()+'.'
	
	m = mrdfits(rawdir+files[pp], 1, head, /silent) ;read in the fits file
	nme = strsplit(files[pp], '.', /extract); get the image name ready
	hextract, m, head, img, nhead, 44, 2091, 0, 2047, /silent;cut out the overscane regions
	sky, img, sky_m0, sky_s0, /silent ;get the initial sky of the image
	sze = size(img);get the image size
	
	;;;;;REMOVE THE FIRST BACKGROUND;;;;

	if grad1 eq 1 then begin
		tmp_bkg = fltarr(sze[1],sze[2], n_elements(bin_ex)) ;make the temporary background holder which iwll be median combined

		;loop over the background to mask stars
		for tt = 0l, n_elements(bin_ex)-1 do begin ;bin loop
			for ii = bin_ex[tt]-bin_box, sze[1]-1, bin_box do begin ; x loop
				for zz = bin_ex[tt]-bin_box, sze[2]-1, bin_box do begin; y loop

					;make sure you dont extend past the image at the end
					mxx = ii+bin_box-1 & mxy = zz+bin_box-1
					if (ii+bin_box-1 gt sze[1]-1) then mxx = sze[1]-1
					if (zz+bin_box-1 gt sze[2]-1) then mxy = sze[2]-1

					;make sure you dont extend past the image in the beginning
					mnx = ii & mny = zz
					if (ii lt 0) then mnx = 0
					if (zz lt 0) then mny = 0

					;select the image in only a box like fashion
					tmp = img[mnx:mxx,mny:mxy]

					;get the sky background in the temporary box
					sky, tmp, sky_m1, sky_s1, /silent

					;only select the pixels which are in a reasonable sky range
					ok1 = where(tmp gt sky_m1-sky_s1 and tmp lt sky_m1+3*sky_s1)

					;assuming there is more than one pixel in a reasonable sky range get the median range
					if n_elements(ok1) gt 1 then begin

						;get the mean clipped background
						meanclip, tmp[ok1], mdn, std, clipsig = 3, maxiter = 100

						;figure out the sky pixels (ok) and the 'stars' (bd)						
						ok = where(tmp lt mdn+3*std, comp = bd)

						;assuming there is at least one 'star' pixel replace it with gaussian sky
						if bd[0] ne -1 then begin
							tmp[bd] = randomn(seed, n_elements(bd))*std+mdn
						endif
			
					endif
					;update the big matrix with the sky background
					tmp_bkg[mnx:mxx,mny:mxy,tt] = tmp
				endfor; end y loop
			endfor; end x loop
		endfor; end bin loop

		;median combine the background images
		medarr, tmp_bkg, fin_bkg

		;smooth the background image with a box x box half the resolution of your initial image
		res = filter_image(fin_bkg, smooth = bin_box/2, /all_pixels)
		print, 'First background level removed at '+systime()+'.'

		;optional write the image if you are debugging etc
		if debug eq 1 then writefits, 'bkg1.fits', res

		;update the image, subtract the sky, but add back the background level
		img_grd = img-res+sky_m0	

		;if debugging write out the cleaned image at this step
		if debug eq 1 then writefits, 'img_grd.fits', img_grd 
		print, 'First background of '+files[pp]+' subtracted at '+systime()+'.'
	endif
	if grad1 eq 0 then img_grd = img

	if strap eq 1 then begin
		;;;;REMOVE THE STRAPS;;;;;

		;get a new background image to remove the straps
		bkg_stp = fltarr(sze[1],sze[2]);make new background image
		sky, img_grd, sky_m2, sky_s2, /silent;get new sky level
		pix = indgen(sze[1]);set up a pixel count for the fitting

		;now loop over each ccd column, and find the straps
		for ii =0l, sze[2]-1 do begin

			;if the flux in a column is larger than a typical flux in a column then remove it
			if (median(img_grd[ii,*]) gt sky_m2+sky_s2) then begin

				res2 = fltarr(sze[2]/bin_stp) ;holder for the background
				pix2 = fltarr(sze[2]/bin_stp) ;holder for the pixels
				uu = 0

				;loop over pixels of size bin_stp to get the median level
				for jj = 0, sze[2]-bin_stp, bin_stp do begin
					mx = jj+bin_stp ;get the max pixel
					if jj+bin_stp gt sze[2]-1 then mx = sze[2]-1 ;make sure the max pixels doesnt extend past image
					resistant_mean, img_grd[ii,jj:mx], 2., m ;get the resistant mean at 2 sigma
					res2[uu] = m & pix2[uu] = jj+bin_stp/2 ;update the hodlers
					uu = uu+1
				endfor
				;interpolate over the pixels to get the strap
				res = interpol(res2, pix2, pix)
				bkg_stp[ii,*] = res;update the background with the straps
			endif
		endfor
		;subtraction the background level from the straps image
		ok = where(bkg_stp ne 0)
		bkg_stp[ok] = bkg_stp[ok]-sky_m2

		;if debugging write the background image
		if debug eq 1 then writefits, 'bkg_stp.fits', bkg_stp 

		;update the image and get new sky level
		img_stp = img_grd-bkg_stp
		sky, img_stp, sky_m3, sky_s3, /silent
	
		;if debugging write the newly cleaned image
		if debug eq 1 then writefits, 'img_stp.fits', img_stp 
		print, 'Straps subtracted from '+files[pp]+' at '+systime()+'.'
	endif
	if strap eq 0 then img_stp = img_grd

	if grad2 eq 1 then begin
		;;; DO ONE LAST GRADIENT BACKGROUND REMOVAL TO CLEAN THINGS UP;;;	
		;make the new background holder
		tmp_bkg2 = fltarr(sze[1],sze[2], n_elements(bin_ex))

		;loop over the background to mask stars
		for tt = 0l, n_elements(bin_ex)-1 do begin ;bin loop
			for ii = bin_ex[tt]-bin_box, sze[1]-1, bin_box do begin ; x loop
				for zz = bin_ex[tt]-bin_box, sze[2]-1, bin_box do begin; y loop

					;make sure you dont extend past the image at the end
					mxx = ii+bin_box-1 & mxy = zz+bin_box-1
					if (ii+bin_box-1 gt sze[1]-1) then mxx = sze[1]-1
					if (zz+bin_box-1 gt sze[2]-1) then mxy = sze[2]-1

					;make sure you dont extend past the image in the beginning
					mnx = ii & mny = zz
					if (ii lt 0) then mnx = 0
					if (zz lt 0) then mny = 0

					;select the image in only a box like fashion
					tmp = img_stp[mnx:mxx,mny:mxy]

					;get the sky background
					sky, tmp, sky_m4, sky_s4, /silent

					;only select the pixels which are in a reasonable sky range
					ok1 = where(tmp gt sky_m4-sky_s4 and tmp lt sky_m4+3*sky_s4)

					;assuming there is more than one pixel in a reasonable sky range get the median range
					if n_elements(ok1) gt 1 then begin
						;get the mean clipped background
						meanclip, tmp[ok1], mdn, std, clipsig = 3, maxiter = 100
						;figure out the sky pixels and the 'stars'						
						ok = where(tmp lt mdn+3*std, comp = bd)
						;assuming there is one 'star' pixel replace it with the sky
						if bd[0] ne -1 then begin
							tmp[bd] = randomn(seed, n_elements(bd))*std+mdn
						endif
				
					endif
					;update the big matrix with the sky background
					tmp_bkg2[mnx:mxx,mny:mxy,tt] = tmp
				endfor
			endfor
		endfor

		;median combine the background images
		medarr, tmp_bkg2, fin_bkg2

		;smooth the backgroudn image with a box x box half the resolution of your initial image
		res2 = filter_image(fin_bkg2, smooth =  bin_box/2, /all_pixels)
		if debug eq 1 then writefits, 'bkg2.fits', res2 ;if debugging is 1 then write out the last gradient image
	
		;subtract the background and add the level back
		fin_img = img_stp-res2+sky_m3
		print, 'Second background of '+files[pp]+' subtracted at '+systime()+'.'
	endif
	if grad2 eq 0 then fin_img = img_stp

	;add stuff to the header
	if grad1 eq 1 then sxaddpar, nhead, 'S_Grad1', 'yes'
	if strap eq 1 then sxaddpar, nhead, 'Strap', 'yes'
	if grad2 eq 1 then sxaddpar, nhead, 'S_Grad2', 'yes'

	;write out the new image
	writefits, clndir+nme[0]+'_s.fits', fin_img, nhead
	print, 'Cleaned '+files[pp]+' written at '+systime()+'.'
endfor

print, 'All done. See ya later alligator.'
END

import main_zf as pr

pr.reload()
pr.zf.fft_gam_2d(pr.dd, {'w_start': 0.5e-1, 'w_end': 1,
                   's_start': 0.5,
                  't_start': 0.5e5})


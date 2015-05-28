path = "C:\\Users\\tew207\\Desktop\\RED_ACCEPTED_FINAL_DATA_CODE\\"

function loadCOVnew96(path::String)

  cov = Array(Float64, (41,41,1))
  Num = similar(cov)

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in77cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in77n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in76cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in76n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in75cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in75n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in74cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in74n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in73cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in73n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in72cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in72n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in71cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in71n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in70cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in70n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh_4yr1in69cov.dta"))
  num = cat(3, num, readdlm(path*"gen96coh_4yr1in69n.dta"))

  cov = cat(3, cov, readdlm(path*"gen96coh1in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh1in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh2in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh2in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh3in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh3in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh4in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh4in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh5in68cov.raw"))
  num = cat(3, num,readdlm(path*"gen96coh5in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh6in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh6in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh7in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh7in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh8in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh8in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh9in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh9in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh10in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh10in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh11in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh11in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh12in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh12in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh13in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh13in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh14in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh14in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh15in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh15in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh16in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh16in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh17in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh17in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh18in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh18in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh19in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh19in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh20in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh20in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh21in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh21in68n.raw"))

  cov = cat(3, cov, readdlm(path*"gen96coh22in68cov.raw"))
  num = cat(3, num, readdlm(path*"gen96coh22in68n.raw"))

  return cov, num
end

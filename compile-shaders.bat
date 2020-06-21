@echo off

REM compiler shaders

REM normal shader
third-party\libs\bgfx\install\bin\shaderc.exe ^
-f shader\next\v_next.sc -o shader\next\v_next.bin ^
--platform windows --type vertex --verbose -i ./ -p vs_5_0

third-party\libs\bgfx\install\bin\shaderc.exe ^
-f shader\next\f_next.sc -o shader\next\f_next.bin ^
--platform windows --type fragment --verbose -i ./ -p ps_5_0

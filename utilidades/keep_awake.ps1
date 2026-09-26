# Impide la suspension del equipo MIENTRAS corra la campana C2.
# No cambia el plan de energia: usa SetThreadExecutionState, que se revierte
# automaticamente al terminar este proceso. Vigila el PID del lanzador y el
# marcador FIN del registro, no el numero de procesos python: entre corrida y
# corrida hay instantes sin ninguno y cortarse ahi seria un falso final.
Add-Type -Name PW -Namespace Win32 -MemberDefinition @'
[DllImport("kernel32.dll")]
public static extern uint SetThreadExecutionState(uint esFlags);
'@
$ES_CONTINUOUS = [uint32]"0x80000000"; $ES_SYSTEM = [uint32]"0x00000001"
$log = "C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\qd_experiments\campana_c2.log"
$lanzador = 18216
while ($true) {
  $fin = (Test-Path $log) -and (Select-String -Path $log -Pattern "FIN:" -Quiet)
  $vivo = @(Get-Process -Id $lanzador -ErrorAction SilentlyContinue).Count -gt 0
  if ($fin -or -not $vivo) { break }
  [Win32.PW]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM) | Out-Null
  Start-Sleep -Seconds 60
}
[Win32.PW]::SetThreadExecutionState($ES_CONTINUOUS) | Out-Null
"keep-awake C2 terminado: " + (Get-Date) |
  Add-Content "C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\qd_experiments\campana_c2.log"

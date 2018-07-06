# echo "cleanup TestAlgo TestAlgo-00-00-01 in /afs/ihep.ac.cn/users/p/privalov/batch/7.0.3"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTestAlgotempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTestAlgotempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt cleanup -csh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  $* >${cmtTestAlgotempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt cleanup -csh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  $* >${cmtTestAlgotempfile}"
  set cmtcleanupstatus=2
  /bin/rm -f ${cmtTestAlgotempfile}
  unset cmtTestAlgotempfile
  exit $cmtcleanupstatus
endif
set cmtcleanupstatus=0
source ${cmtTestAlgotempfile}
if ( $status != 0 ) then
  set cmtcleanupstatus=2
endif
/bin/rm -f ${cmtTestAlgotempfile}
unset cmtTestAlgotempfile
exit $cmtcleanupstatus


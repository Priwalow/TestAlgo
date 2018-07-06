# echo "setup TestAlgo TestAlgo-00-00-01 in /afs/ihep.ac.cn/users/p/privalov/batch/7.0.3"

if ( $?CMTROOT == 0 ) then
  setenv CMTROOT /afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/contrib/CMT/v1r25
endif
source ${CMTROOT}/mgr/setup.csh
set cmtTestAlgotempfile=`${CMTROOT}/mgr/cmt -quiet build temporary_name`
if $status != 0 then
  set cmtTestAlgotempfile=/tmp/cmt.$$
endif
${CMTROOT}/mgr/cmt setup -csh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  -no_cleanup $* >${cmtTestAlgotempfile}
if ( $status != 0 ) then
  echo "${CMTROOT}/mgr/cmt setup -csh -pack=TestAlgo -version=TestAlgo-00-00-01 -path=/afs/ihep.ac.cn/users/p/privalov/batch/7.0.3  -no_cleanup $* >${cmtTestAlgotempfile}"
  set cmtsetupstatus=2
  /bin/rm -f ${cmtTestAlgotempfile}
  unset cmtTestAlgotempfile
  exit $cmtsetupstatus
endif
set cmtsetupstatus=0
source ${cmtTestAlgotempfile}
if ( $status != 0 ) then
  set cmtsetupstatus=2
endif
/bin/rm -f ${cmtTestAlgotempfile}
unset cmtTestAlgotempfile
exit $cmtsetupstatus


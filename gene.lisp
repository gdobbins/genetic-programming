(in-package :cl-user)
(eval-when (:compile-toplevel :load-toplevel :execute)
  ;(ql:quickload :series)
  (ql:quickload :alexandria))
(defpackage :genetic-programming
  (:use :common-lisp :alexandria)
  (:export :gene-prog :make-population :population-genes :population-gene-fitness-list :population-gene-adjusted-fitness-list :population-number-of-genes :population-average-fitness :population-maximum-fitness :population-minimum-fitness :population-best-gene-so-far :population-best-gene-so-far-fitness :make-functions :functions-function-list :functions-arg-num-list))
(in-package :genetic-programming)

;(defclass gene ()
;  ((gene :accessor gene :initarg :gene)
;   (raw-fitness :accessor raw-fitness)
;   (adjusted-fitness :accessor adjusted-fitness)))

(defstruct population
  (genes nil :type (or cons null))
  (gene-fitness-list nil :type (or cons null))
  (gene-adjusted-fitness-list nil :type (or cons null))
  (number-of-genes 0 :type fixnum)
  (average-fitness 0 :type number)
  (maximum-fitness 0 :type number)
  (minimum-fitness 0 :type number)
  (best-gene-so-far nil :type (or cons null))
  (best-gene-so-far-fitness 0 :type number))

(defstruct functions
  (function-list nil :type (or cons null))
  (arg-num-list nil :type (or cons null)))

(defun tree-replace (tree find max-depth min-depth replace &rest vars)
  (loop for i in tree collect
       (if (and (consp i) (<= 0 max-depth) (not #1=(find-if (lambda (x) (equalp i x)) (last vars))))
	   (apply 'tree-replace i find (1- max-depth) (1- min-depth) replace vars)
	   (if (or (equal i find) (and (<= 0 min-depth) #1#) (and (>= 0 max-depth) (consp i) (not #1#)))
	       (apply replace vars)
	       i))))

(defun tree-find-not (tree find)
  (let ((rep t))
    (loop for i in tree do
	 (if rep
	     (if (consp i)
		 (setf rep (tree-find-not i find))
		 (if (equal i find)
		     (setf rep nil)))))
    rep))

(defgeneric create-random-gene (functions terminals &optional max-depth min-depth))

(defmethod create-random-gene ((functions functions) (terminals cons) &optional (max-depth 6) (min-depth 2))
  (let* ((len-fun (length #1=(functions-function-list functions))) (len-term (length terminals)) (len (+ len-fun len-term)) gene)
    (flet ((create-new-sexp (num functions terminals)
		 (if (< num len-fun)
		     (append (list (nth num #1#)) (make-list (nth num (functions-arg-num-list functions)) :initial-element :sexp))
		     (nth (- num len-fun) terminals))))
      (setf gene (create-new-sexp (random len-fun) functions terminals))
      (do () ((tree-find-not gene :sexp) gene)
	(setf gene (tree-replace gene :sexp max-depth min-depth (lambda (x y z) (create-new-sexp (random x) y z)) len functions terminals))))))

(defun initialize-population (current-population functions terminals generation)
  (declare (type population current-population) (type cons terminals))
  (assert (= 0 generation))
  (if (< (length #1=(population-genes current-population)) #2=(population-number-of-genes current-population))
      (tagbody (setf (population-genes current-population) (append (list (list (create-random-gene functions terminals))) (population-genes current-population))) (initialize-population current-population functions terminals generation))
;      (initialize-population (make-population :genes (append (list (create-random-gene functions terminals)) #1#) :number-of-genes #2#) functions terminals generation)
      current-population))

(defun adjust-to-linear-fitness (fitness-list max-fit min-fit avg-fit &optional (max-fit-multiplier 2))
  (declare (type cons fitness-list) (type number max-fit min-fit avg-fit max-fit-multiplier))
  (if (= max-fit avg-fit) (make-list (length fitness-list) :initial-element 1)
      (let ((a 0) (b 0) result)
	(setf a (/ (- max-fit-multiplier 1) (- max-fit avg-fit)) b (* -1 (1+ (* a avg-fit))))
	#1=(setf result (loop for i in fitness-list collect (+ b (* a i))))
	(if (some (lambda (x) (> 0 x)) result)
	    (if (= avg-fit min-fit)
		(setf result (make-list (length fitness-list) :initial-element 1))
		(and (setf a (/ 1 (- avg-fit min-fit)) b (* -1 a min-fit)) #1#)))
	result)))

(defun update-population-fitness (current-population fitness-function)
  (declare (type function fitness-function))
  (labels ((upf-helper (gene-list fitness-list current-population fitness-function)
	     (declare (type function fitness-function))
	     (if gene-list
		 (tagbody (setf (first fitness-list) (funcall fitness-function (first gene-list))) (upf-helper (rest gene-list) (rest fitness-list) current-population fitness-function))))
	   (upf-best-gene-helper (gene-list fitness-list best-gene-fitness)
	     (if (= (car fitness-list) best-gene-fitness)
		 (values (car gene-list) (car fitness-list))
		 (upf-best-gene-helper (cdr gene-list) (cdr fitness-list) best-gene-fitness))))
    (upf-helper (population-genes current-population) #1=(population-gene-fitness-list current-population) current-population fitness-function)
    (setf (population-average-fitness current-population) (/ (loop for i in #1# summing i) (population-number-of-genes current-population)))
    (setf (population-minimum-fitness current-population) (loop for i in #1# minimize i))
    (setf (population-maximum-fitness current-population) (loop for i in #1# maximize i))
    (setf (population-gene-adjusted-fitness-list current-population) (adjust-to-linear-fitness (population-gene-fitness-list current-population) (population-maximum-fitness current-population) (population-minimum-fitness current-population) (population-average-fitness current-population)))
    (if (> (population-maximum-fitness current-population) (population-best-gene-so-far-fitness current-population))
	(multiple-value-bind (gene fitness) (upf-best-gene-helper (population-genes current-population) (population-gene-fitness-list current-population) (population-maximum-fitness current-population))
	  (setf (population-best-gene-so-far current-population) gene (population-best-gene-so-far-fitness current-population) fitness)))))

(defun create-mating-pool (current-population)
  (let* ((total-fitness (coerce (loop for i in (population-gene-adjusted-fitness-list current-population) sum i) 'float)) (num-genes (population-number-of-genes current-population)))
    (labels ((cmp-helper (gene-list fitness-list selector-number)
	       (if (> (car fitness-list) selector-number)
		   (car gene-list)
		   (cmp-helper (cdr gene-list) (cdr fitness-list) (- selector-number (car fitness-list))))))
      (make-population :number-of-genes num-genes :best-gene-so-far (population-best-gene-so-far current-population) :best-gene-so-far-fitness (population-best-gene-so-far-fitness current-population)
		       :average-fitness (population-average-fitness current-population)
		       :maximum-fitness (population-maximum-fitness current-population)
		       :minimum-fitness (population-minimum-fitness current-population)
		       :gene-fitness-list (population-gene-fitness-list current-population)
		       :genes (loop for i from 1 to num-genes collect (copy-tree (cmp-helper (population-genes current-population) (population-gene-adjusted-fitness-list current-population) (random total-fitness))))))))

(defun breed-mating-pool (mating-pool functions terminals &optional (breed-percentage 0.9) (mutate-percentage .001))
  (let ((num-breeding-pairs (round (* (population-number-of-genes mating-pool) breed-percentage))))
    (labels ((bmp-helper (gene-list functions terminals pairs-remaining mutate-percentage)
	       (if (<= 2 pairs-remaining)
		   (tagbody (rotatef (car (random-elt (node-list (first gene-list)))) (car (random-elt (node-list (second gene-list)))))
		      (bmp-helper (cddr gene-list) functions terminals (- pairs-remaining 2) mutate-percentage))
		   (let ((n (length gene-list)))
		     (dotimes (i (round (+ (* n mutate-percentage (- 1 mutate-percentage) (gaussian-random)) (* n mutate-percentage))))
		       (setf (car (random-elt (node-list (random-elt gene-list)))) (create-random-gene functions terminals 3 1)))))))
      (bmp-helper (population-genes mating-pool) functions terminals num-breeding-pairs mutate-percentage))
    mating-pool))

(defun update-history-report (current-population generation print-freq)
  (if (zerop (mod generation print-freq))
      (format t "~&--------------------------------------------------~%Generation: ~d     Stardate: ~{~d-~d-~d  ~d:~2,'0d:~2,'0d~}~%Maximum Fitness: ~d   Minimum Fitness: ~d   Average Fitness: ~d   Best Fitness: ~d   "
	      generation (nreverse (subseq (multiple-value-list (get-decoded-time)) 0 6)) (population-maximum-fitness current-population) (population-minimum-fitness current-population) (population-average-fitness current-population) (population-best-gene-so-far-fitness current-population))))

(let ((current-population (make-population)) (mating-pool (make-population))
      (terminals nil) (functions nil) (print-freq 1)
      (fitness-function #'(lambda (&rest rest) (declare (ignore rest)) 1))
      (generation 0))
  (defun gene-prog (&key (operation :next-generation) (input) (generations 50);for progressing to the next evolution
		      (delete-duplicates nil) (population-size 100) (init-functions functions) (init-terminals terminals) (fit-function fitness-function) (print-frequency print-freq);for during initialization, whether to check for duplicates
		      (file nil)) ;for saving and loading
    (declare (type fixnum generation))
    (case operation
      (:next-generation
       (update-population-fitness current-population fit-function)
       (setf mating-pool (create-mating-pool current-population))
       (setf current-population (breed-mating-pool mating-pool functions terminals))
       (update-history-report current-population generation print-freq)
       (setf generation (1+ generation)))
      (:best-of-generation
       (labels ((bog-helper (gene-list fitness-list max-fitness)
		  (if (= max-fitness (car fitness-list))
		      (car gene-list)
		      (bog-helper (cdr gene-list) (cdr fitness-list) max-fitness))))
	 (update-population-fitness current-population fit-function)
	 (values (bog-helper (population-genes current-population) (population-genes current-population) (population-maximum-fitness current-population)) (population-maximum-fitness current-population))))
      (:best-of-generation-continue
       (setf mating-pool (create-mating-pool current-population))
       (setf current-population (breed-mating-pool mating-pool functions terminals))
       (update-history-report current-population generation print-freq)
       (setf generation (1+ generation)))
      (:next-generations
       (setf fitness-function fit-function)
       (dotimes (i generations)
	 (update-population-fitness current-population fitness-function)
	 (setf mating-pool (create-mating-pool current-population))
	 (setf current-population (breed-mating-pool mating-pool functions terminals))
	 (update-history-report current-population generation print-freq)
	 (setf generation (1+ generation))))
      (:change-terminals (setf terminals init-terminals))
      (:initialize
       (setf terminals init-terminals functions init-functions
	     (population-number-of-genes current-population) population-size
	     (population-number-of-genes mating-pool) population-size)
       (setf fitness-function fit-function)
       (setf print-freq print-frequency)
       (initialize-population current-population functions terminals generation)
       (setf (population-gene-fitness-list current-population) (make-list (population-number-of-genes current-population) :initial-element 0))
       (if delete-duplicates (funcall 'gene-prog :operation :delete-initial-duplicates
				      :population-size population-size
				      :init-functions init-functions
				      :init-terminals init-terminals
				      :delete-duplicates t)))
      (:delete-initial-duplicates
       (let ((old-population-size (length (population-genes current-population))))
       (setf #1=(population-genes current-population) (delete-duplicates #1# :test #'tree-equal))
       (funcall 'gene-prog :operation :initialize
		:population-size old-population-size
		:init-functions init-functions
		:init-terminals init-terminals
		:print-frequency print-frequency
		:delete-duplicates (if (= old-population-size (length #1#)) nil t))))
      (:save (with-open-file (out file :direction :output :if-exists :overwrite)
	       (prin1 current-population out)
	       (prin1 mating-pool out)
	       (prin1 terminals out)
	       (prin1 functions out)
	       (prin1 generation out)
	       (prin1 print-freq out)))
      (:load (with-open-file (in file :direction :input)
	       (setf current-population (read in))
	       (setf mating-pool (read in))
	       (setf terminals (read in))
	       (setf functions (read in))
	       (setf generation (read in))
	       (setf print-freq (read in))))
      (:set-fit-function (setf fitness-function fit-function))
      (:mating-pool mating-pool)
      (:fitness-function fitness-function)
      (:input (values input generations))
      (:terminals terminals)
      (:functions functions)
      (:generation generation)
      (:print-freq print-freq)
      (:current-population current-population))))

(defun node-list (gene &optional (acc nil) (ret nil))
  (declare (type (or cons null) gene acc ret))
  (if gene
      (if (consp (car gene))
	  (node-list (cdar gene) (cons gene acc) (cons (cdr gene) ret))
	  (node-list (cdr gene) (cons gene acc) ret))
      (if ret
	  (node-list (first ret) acc (rest ret))
	  acc)))

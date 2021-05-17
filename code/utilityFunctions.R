# cosineSim computes the standard cosine similarity of two vectors: https://en.wikipedia.org/wiki/Cosine_similarity
# the vectors can be of any length; they do not have to share indices.
# each vector has a set of values, i.e., vec1 and vec2 do not contain duplicate values. 
# For example, 
#cosineSim(c(1,3,5), c(4,9,7,6)) = 0.0
#cosineSim(c(1,3,5),c(1,3,5)) = 1.0
#cosineSim(c(1,5),c(1,3,5)) = 0.81
#cosineSim(c(1,7,9,5),c(1,3,5))=0.57
cosineSim <- function(vec1, vec2) {
   len1 = length(vec1)
   len2 = length(vec2)
   result = 0.0
   if (is.na(vec1) || len1 == 0 ||
       is.na(vec2) || len2 == 0) {
      return(0.0)
   }
   
   inter = length(intersect(vec1, vec2))
   if (inter == 0) {
      result = 0.0
   }
   else{
      result = inter / (sqrt(len1) * sqrt(len2))
   }
   
   return(result)
}
